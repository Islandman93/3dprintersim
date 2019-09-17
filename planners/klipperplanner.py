# Code for simulating klipper planner, most of it is copied from klipper/klippy/toolhead.py
# Original license below
#
# Copyright (C) 2016-2018  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import math
import numpy as np
from .util import find_max_vel
LOOKAHEAD_FLUSH_TIME = 0.250
DA_PERIOD = 0.032


class KlipperPlanner:
    def __init__(self, acceleration, max_vel, square_corner_velocity, accel_to_decel=None, dynamic_acceleration=False):
        self.acceleration = acceleration
        self.da_min_accel = acceleration / 2
        self.da_min_period = DA_PERIOD*2
        self.do_dynamic_accel = dynamic_acceleration
        if accel_to_decel is None:
            self.accel_to_decel = acceleration / 2
        else:
            self.accel_to_decel = accel_to_decel
        # assumes max_vel is total equating to how much plastic can be extruded
        self.max_vel = max_vel
        scv2 = square_corner_velocity**2
        self.junction_deviation = scv2 * (math.sqrt(2.) - 1.) / self.acceleration

        # start copy from klipper
        self.queue = []
    def flush(self):
        queue = self.queue
        flush_count = len(queue)
        update_flush_count = False
        # Traverse queue from last to first move and determine maximum
        # junction speed assuming the robot comes to a complete stop
        # after the last move.
        delayed = []
        next_end_v2 = next_smoothed_v2 = peak_cruise_v2 = 0.
        for i in range(flush_count-1, 0, -1):
            move = queue[i]
            reachable_start_v2 = next_end_v2 + move.delta_v2
            start_v2 = min(move.max_start_v2, reachable_start_v2)
            reachable_smoothed_v2 = next_smoothed_v2 + move.smooth_delta_v2
            smoothed_v2 = min(move.max_smoothed_v2, reachable_smoothed_v2)
            if smoothed_v2 < reachable_smoothed_v2:
                # It's possible for this move to accelerate
                if (smoothed_v2 + move.smooth_delta_v2 > next_smoothed_v2
                    or delayed):
                    # This move can decelerate or this is a full accel
                    # move after a full decel move
                    if update_flush_count and peak_cruise_v2:
                        flush_count = i
                        update_flush_count = False
                    peak_cruise_v2 = min(move.max_cruise_v2, (
                        smoothed_v2 + reachable_smoothed_v2) * .5)
                    if delayed:
                        # Propagate peak_cruise_v2 to any delayed moves
                        if not update_flush_count and i < flush_count:
                            for m, ms_v2, me_v2 in delayed:
                                mc_v2 = min(peak_cruise_v2, ms_v2)
                                m.set_junction(min(ms_v2, mc_v2), mc_v2
                                               , min(me_v2, mc_v2))
                        del delayed[:]
                if not update_flush_count and i < flush_count:
                    cruise_v2 = min((start_v2 + reachable_start_v2) * .5
                                    , move.max_cruise_v2, peak_cruise_v2)
                    move.set_junction(min(start_v2, cruise_v2), cruise_v2
                                      , min(next_end_v2, cruise_v2))
            else:
                # Delay calculating this move until peak_cruise_v2 is known
                delayed.append((move, start_v2, next_end_v2))
            next_end_v2 = start_v2
            next_smoothed_v2 = smoothed_v2
        if update_flush_count:
            return
        # Generate step times for all moves ready to be flushed
        # Simulation: first move is not set, do that here
        move_stuff = []
        # TODO: sometimes second move doesn't have start v set
        end_vel_first_move = queue[1].start_v
        max_vel_first_move = find_max_vel(queue[0].move_d, 0, end_vel_first_move, self.accel_to_decel)
        cruise_vel2_first_move = min(max_vel_first_move**2, self.max_vel**2)
        queue[0].set_junction(0, cruise_vel2_first_move, end_vel_first_move**2)
        # Remove processed moves from the queue
        return queue
    def add_move(self, start_pos, end_pos):
        move = self.generateMove(start_pos, end_pos)
        self.queue.append(move)
        if len(self.queue) == 1:
            return
        move.calc_junction(self.queue[-2])

    def generateMove(self, start_pos, end_pos):
        if self.do_dynamic_accel:
            minimum_da_move_d = 0.5 * self.da_min_accel * self.da_min_period**2
        else:
            minimum_da_move_d = float('inf')

        # calc distance
        axes_d = [end_pos[i] - start_pos[i] for i in (0, 1)] + [0]
        move_d = math.sqrt(sum([d*d for d in axes_d[:2]]))

        # if less than minimum dynamic acceleration use regular move
        if move_d < minimum_da_move_d:
            move = Move(start_pos, end_pos, self.max_vel, self.acceleration, self.accel_to_decel, self.junction_deviation)
        # else possible to use dynamic acceleration
        else:
            move = DynamicAccelerationMove(start_pos, end_pos, self.max_vel, self.acceleration, self.accel_to_decel, self.junction_deviation)
        return move


class Move:
    def __init__(self, start_pos, end_pos, speed, acceleration, accel_to_decel, junction_deviation):
        self.start_pos = tuple(start_pos)
        self.end_pos = tuple(end_pos)
        self.accel = acceleration
        self.decel = acceleration
        self.accel_to_decel = accel_to_decel
        self.junction_deviation = junction_deviation
        self.junction_accel = acceleration
        velocity = speed
        self.axes_d = axes_d = [end_pos[i] - start_pos[i] for i in (0, 1)] + [0]
        self.move_d = move_d = math.sqrt(sum([d*d for d in axes_d[:2]]))
        self.min_move_t = move_d / velocity
        # Junction speeds are tracked in velocity squared.  The
        # delta_v2 is the maximum amount of this squared-velocity that
        # can change in this move.
        self.max_start_v2 = 0.
        self.max_cruise_v2 = velocity**2
        self.delta_v2 = 2.0 * move_d * self.accel
        self.max_smoothed_v2 = 0.
        self.smooth_delta_v2 = 2.0 * move_d * self.accel_to_decel
    def calc_junction(self, prev_move):
        # Find max velocity using "approximated centripetal velocity"
        axes_d = self.axes_d
        prev_axes_d = prev_move.axes_d
        junction_cos_theta = -((axes_d[0] * prev_axes_d[0]
                                + axes_d[1] * prev_axes_d[1]
                                + axes_d[2] * prev_axes_d[2])
                               / (self.move_d * prev_move.move_d))
        if junction_cos_theta > 0.999999:
            return
        junction_cos_theta = max(junction_cos_theta, -0.999999)
        sin_theta_d2 = math.sqrt(0.5*(1.0-junction_cos_theta))
        R = (self.junction_deviation * sin_theta_d2
             / (1. - sin_theta_d2))
        tan_theta_d2 = sin_theta_d2 / math.sqrt(0.5*(1.0+junction_cos_theta))
        move_centripetal_v2 = .5 * self.move_d * tan_theta_d2 * self.junction_accel
        prev_move_centripetal_v2 = (.5 * prev_move.move_d * tan_theta_d2
                                    * prev_move.junction_accel)
        self.max_start_v2 = min(
            R * self.accel, R * prev_move.accel,
            move_centripetal_v2, prev_move_centripetal_v2,
            self.max_cruise_v2, prev_move.max_cruise_v2,
            prev_move.max_start_v2 + prev_move.delta_v2)
        self.max_smoothed_v2 = min(
            self.max_start_v2
            , prev_move.max_smoothed_v2 + prev_move.smooth_delta_v2)

    def set_junction(self, start_v2, cruise_v2, end_v2):
        # Determine accel, cruise, and decel portions of the move distance
        inv_delta_v2 = 1. / self.delta_v2
        self.accel_r = accel_r = (cruise_v2 - start_v2) * inv_delta_v2
        self.decel_r = decel_r = (cruise_v2 - end_v2) * inv_delta_v2
        self.cruise_r = cruise_r = 1. - accel_r - decel_r
        # Determine move velocities
        self.start_v = start_v = math.sqrt(start_v2)
        self.cruise_v = cruise_v = math.sqrt(cruise_v2)
        self.end_v = end_v = math.sqrt(end_v2)
        # Determine time spent in each portion of move (time is the
        # distance divided by average velocity)
        self.accel_t = accel_r * self.move_d / ((start_v + cruise_v) * 0.5)
        self.cruise_t = cruise_r * self.move_d / cruise_v
        self.decel_t = decel_r * self.move_d / ((end_v + cruise_v) * 0.5)

        # if dynamic acceleration
        da_accel, da_accel_t = None, None
        if cruise_v - start_v > 0:
            # TODO: min dynamic accel
            da_accel = max((cruise_v-start_v) / DA_PERIOD, 1)
            da_accel_t = (cruise_v-start_v) / da_accel
        da_decel, da_decel_t = None, None
        if cruise_v - end_v > 0:
            # TODO: min dynamic accel
            da_decel = max((cruise_v-end_v) / DA_PERIOD, 1)
            da_decel_t = (cruise_v-end_v) / da_decel

        # if self.accel_t > DA_PERIOD:
            # print('missed accel')
            # self.accel = (cruise_v - start_v) / DA_PERIOD
            # self.accel_t = DA_PERIOD
            # changed = True
        # if self.decel_t > DA_PERIOD:
            # print('missed decel')
            # self.decel = (cruise_v - end_v) / DA_PERIOD
            # self.decel_t = DA_PERIOD
            # changed = True
        # recalculate times
        # if da_accel is not None and da_decel is not None:
        if False:
            # Determine time spent in cruise
            accel_dist = .5 * (start_v + cruise_v) * da_accel_t
            decel_dist = .5 * (cruise_v + end_v) * da_decel_t
            cruise_dist = self.move_d - (accel_dist + decel_dist)
            # da works
            if cruise_dist > 0:
                self.accel = da_accel
                self.accel_t = da_accel_t
                self.decel = da_decel
                self.decel_t = da_decel_t
                self.cruise_t = cruise_dist / cruise_v
            # else:
                # print('revert')

    def __str__(self):
        return "{:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(self.start_v, self.cruise_v, self.end_v, self.accel_t, self.cruise_t, self.decel_t, self.delta_v2)


class DynamicAccelerationMove(Move):
    def __init__(self, start_pos, end_pos, speed, acceleration, accel_to_decel, junction_deviation):
        super().__init__(start_pos, end_pos, speed, acceleration, accel_to_decel, junction_deviation)
        self.da_period = DA_PERIOD
        self.da_min_period = self.da_period * 2
        self.da_min_accel = acceleration / 2
        # Junction speeds are tracked in velocity squared.  The
        # delta_v2 is the maximum amount of this squared-velocity that
        # can change in this move.
        self.junction_accel = self.da_min_accel
        self.delta_v2 = 2.0 * self.move_d * self.da_min_accel
        self.smooth_delta_v2 = self.delta_v2

    def set_junction(self, start_v2, cruise_v2, end_v2):
        # Determine move velocities
        self.start_v = start_v = math.sqrt(start_v2)
        self.cruise_v = cruise_v = math.sqrt(cruise_v2)
        self.end_v = end_v = math.sqrt(end_v2)
        # Determine accel, cruise, and decel portions of the move distance
        delta_cruise_start = cruise_v - start_v
        if delta_cruise_start > 0:
            self.accel = max(delta_cruise_start / self.da_period, self.da_min_accel)
            self.accel_t = delta_cruise_start / self.accel
        else:
            self.accel = 0
            self.accel_t = 0
        delta_cruise_end = abs(cruise_v - end_v)
        if delta_cruise_end > 0:
            self.decel = max(delta_cruise_end / self.da_period, self.da_min_accel)
            self.decel_t = delta_cruise_end / self.decel
        else:
            self.decel = 0
            self.decel_t = 0
        # Determine time spent in cruise
        accel_dist = .5 * (start_v + cruise_v) * self.accel_t
        decel_dist = .5 * (cruise_v + end_v) * self.decel_t
        cruise_dist = self.move_d - (accel_dist + decel_dist)
        self.cruise_t = cruise_dist / cruise_v
        if cruise_dist < -0.0001:
            print('cruise error')
            print('req dist', self.move_d, accel_dist, decel_dist)
            print('req vel', self.start_v, self.cruise_v, self.end_v)
            print('calc v', math.sqrt(self.delta_v2))
            print('calc accel', self.accel, self.decel)
            print('calc times', self.accel_t, self.decel_t)

