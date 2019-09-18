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
    def __init__(self, acceleration, max_vel, square_corner_velocity, accel_to_decel=None, acceleration_type='constant', jerk=0):
        self.acceleration = acceleration
        self.da_min_accel = acceleration / 2
        self.da_min_period = DA_PERIOD*2
        self.acceleration_type = acceleration_type
        # only used when jerk limited moves are used
        self.jerk = jerk
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
            reachable_start_v2 = move.calc_max_v2(next_end_v2, move.accel)
            start_v2 = min(move.max_start_v2, reachable_start_v2)
            reachable_smoothed_v2 = move.calc_max_v2(next_smoothed_v2
                    , move.smoothed_accel)
            smoothed_v2 = min(move.max_smoothed_v2, reachable_smoothed_v2)
            if smoothed_v2 < reachable_smoothed_v2:
                # It's possible for this move to accelerate
                if (move.calc_max_v2(smoothed_v2, move.smoothed_accel)
                        > next_smoothed_v2 or delayed):
                    # This move can decelerate or this is a full accel
                    # move after a full decel move
                    if update_flush_count and peak_cruise_v2:
                        flush_count = i
                        update_flush_count = False
                    peak_cruise_v2 = min(move.max_cruise_v2,
                                         move.calc_peak_v2(smoothed_v2,
                                                           next_smoothed_v2,
                                                           move.smoothed_accel))
                    if delayed:
                        # Propagate peak_cruise_v2 to any delayed moves
                        if not update_flush_count and i < flush_count:
                            for m, ms_v2, me_v2 in delayed:
                                mc_v2 = min(peak_cruise_v2, ms_v2)
                                m.set_junction(min(ms_v2, mc_v2), mc_v2
                                               , min(me_v2, mc_v2))
                        del delayed[:]
                if not update_flush_count and i < flush_count:
                    cruise_v2 = min(move.max_cruise_v2, peak_cruise_v2,
                                    move.calc_peak_v2(
                                        start_v2, next_end_v2, move.accel)
                                    )
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
        if self.acceleration_type == 'dynamic' or self.acceleration_type == 'destructive':
            minimum_da_move_d = 0.5 * self.da_min_accel * self.da_min_period**2
        else:
            minimum_da_move_d = float('inf')

        # calc distance
        axes_d = [end_pos[i] - start_pos[i] for i in (0, 1)] + [0]
        move_d = math.sqrt(sum([d*d for d in axes_d[:2]]))

        # if less than minimum dynamic acceleration use regular move
        if move_d < minimum_da_move_d:
            if 'jerk-limited' in self.acceleration_type:
                move = JerkLimitedSmoothMove(start_pos, end_pos, self.max_vel,
                                             self.acceleration, self.accel_to_decel,
                                             self.junction_deviation, jerk=self.jerk)
            else:
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
        self.smoothed_accel = accel_to_decel
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

    def calc_peak_v2(self, start_v2, end_v2, accel):
        return (start_v2 + end_v2 + 2 * self.move_d * accel) * 0.5

    def calc_max_v2(self, start_v2, accel, dist=None):
        dist = dist or self.move_d
        # Check if accel is the limiting factor
        return start_v2 + 2.0 * dist * accel

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


class JerkLimitedSmoothMove(Move):
    def __init__(self, *args, jerk=10):
        super().__init__(*args)
        self.jerk = jerk

    def calc_junction(self, prev_move):
        self.prev_move = prev_move
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
        move_centripetal_v2 = .5 * self.move_d * tan_theta_d2 * self.accel
        prev_move_centripetal_v2 = (.5 * prev_move.move_d * tan_theta_d2
                                    * prev_move.accel)
        self.max_start_v2 = min(
            R * self.accel, R * prev_move.accel,
            move_centripetal_v2, prev_move_centripetal_v2,
            self.max_cruise_v2, prev_move.max_cruise_v2,
            prev_move.calc_max_v2(prev_move.max_start_v2, prev_move.accel))
        self.max_smoothed_v2 = min(
            self.max_start_v2,
            prev_move.calc_max_v2(prev_move.max_smoothed_v2
                , prev_move.smoothed_accel))

    def calc_max_v2(self, start_v2, accel, dist=None):
        dist = dist or self.move_d
        # Check if accel is the limiting factor
        max_accel_v2 = start_v2 + 2.0 * dist * accel
        # Compute maximum achievable speed with limited kinematic jerk using
        # max(jerk) == 6 * accel / accel_time, which is exact for accel order 4
        # and is quite accurate for accel order 6:
        # max(jerk) == 10 / sqrt(3) * accel / accel_time ~=
        #     5.774 * accel / accel_time
        # This leads to the cubic equation
        # (max_v^2 - start_v^2) * (max_v + start_v) / 2 ==
        #     dist^2 * jerk / 3
        # which is solved using Cardano's formula.
        start_v = math.sqrt(start_v2)
        a = 2./3. * start_v
        b = a*a*a
        c = dist * dist * self.jerk / 3.
        if b * 54 < c:
            # Make max_v monotonic over start_v: return the max velocity
            # which works for any start_v velocity below the threshold.
            max_v = 1.5 * (c*.5)**(1./3.)
        else:
            d = math.sqrt(c * (c + 2 * b))
            e = (b + c + d)**(1./3.)
            if e < 0.000000001:
                return start_v
            max_v = e + a*a / e - start_v / 3.
        return min(max_v * max_v, max_accel_v2)

    def calc_effective_accel(self, start_v, cruise_v):
        effective_accel = min(
                self.accel, math.sqrt(self.jerk * (cruise_v - start_v) / 6.))
        return effective_accel

    def calc_min_accel_time(self, start_v, cruise_v):
        min_accel_time = (cruise_v - start_v) / self.accel
        min_accel_time = max(min_accel_time
                , math.sqrt(6. * (cruise_v - start_v) / self.jerk))
        return min_accel_time

    def calc_peak_v2(self, start_v2, end_v2, accel):
        peak_v2 = (start_v2 + end_v2 + 2 * self.move_d * accel) * 0.5
        start_v = math.sqrt(start_v2)
        end_v = math.sqrt(end_v2)
        # This is only an approximate of the min distance for acceleration.
        min_accel_d = self.calc_min_accel_time(min(start_v, end_v)
                , max(start_v, end_v)) * (start_v + end_v) * 0.5
        extra_d = (self.move_d - min_accel_d) * 0.5
        peak_v2 = min(peak_v2
                , self.calc_max_v2(max(start_v2, end_v2), accel, dist=extra_d)
                # With accel_order > 2 this term can get smaller than
                # max(start_v2, end_v2) if min_accel_d is insufficient.
                , self.calc_max_v2(min(start_v2, end_v2), accel
                    , dist=self.move_d-extra_d))
        return max(start_v2, end_v2, peak_v2)

    def set_junction(self, start_v2, cruise_v2, end_v2):
        # Determine move velocities
        self.start_v = start_v = math.sqrt(start_v2)
        self.cruise_v = cruise_v = math.sqrt(cruise_v2)
        self.end_v = end_v = math.sqrt(end_v2)
        # Determine the effective accel and decel
        self.effective_accel = self.calc_effective_accel(start_v, cruise_v)
        self.effective_decel = self.calc_effective_accel(end_v, cruise_v)
        # Determine time spent in each portion of move (time is the
        # distance divided by average velocity)
        self.accel_t = accel_t = self.calc_min_accel_time(start_v, cruise_v)
        self.decel_t = decel_t = self.calc_min_accel_time(end_v, cruise_v)
        self.cruise_t = cruise_t = (self.move_d
                - accel_t * (start_v + cruise_v) * 0.5
                - decel_t * (end_v + cruise_v) * 0.5) / cruise_v
        self.accel = self.effective_accel
        self.decel = self.effective_decel
        if cruise_t < -0.000000001:
            raise Exception(
                    'Logic error: impossible move ms_v2=%.3lf, mc_v2=%.3lf'
                    ', me_v2=%.3lf with move_d=%.3lf, accel=%.3lf, jerk=%.3lf'
                    % (start_v2, cruise_v2, end_v2, self.move_d, self.accel
                        , self.jerk))

