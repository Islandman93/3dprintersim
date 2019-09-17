import math
from planners.util import decompose_xy, smoothstep_accel, smootherstep_accel, check_dynamic_acceleration, stiffness_from_period_damping, destruct

import pymunk as pm
from pymunk import Vec2d


class Simulation:
    def __init__(self, num_axis, periods, damping_coefs=(0.1,0.1), sample_rate=40000, eps=0.0125):
        # Physics setup
        self.space = pm.Space()
        self.space.damping = 0.99
        self.sample_rate = sample_rate
        self.COLLTYPE_DEFAULT = 0
        self.rest_length = 100

        # create axes and springs
        self.axes = []
        self.periods = periods
        self.stiffnesses = []
        pos = 0
        for a in range(num_axis):
            stiffness = stiffness_from_period_damping(periods[a], damping_coefs[a])
            self.stiffnesses.append(stiffness)
            self.axes.append(self.create_axis(pos, stiffness, damping_coefs[a]))
            pos += 100

        # metrics for tracking error on each axis
        self.axes_metrics = []
        for i in range(len(self.axes)):
            self.axes_metrics.append({
                'accel': [],
                'expected': [],
                'actual': [],
                'error': []
            })

        # in some cases float precision on moves along with jerk/junction deviation algs
        # create an impossible move where the desired distance is smaller than the commanded
        # just use some epsilon for really small distances (in mm)
        self.eps = eps

        # for dynamic acceleration, counts the number of unsuccessful applications
        self.da_revert = 0

    def move_xy(self, xy_dist, acceleration, start_vel, max_vel, end_vel, accel_type):
        # save start position for accuracy check at the end of this function
        start_pos = [a[0].position[0] for a in self.axes]
        # assumes max_vel is total equating to how much plastic can be extruded
        # decompose x/y amounts of move, accel, vels
        x_ratio, y_ratio = decompose_xy(xy_dist[0], xy_dist[1])
        deceleration = acceleration
        accel_type, decel_type = accel_type, accel_type

        # dist
        dist = math.sqrt(xy_dist[0] ** 2 + xy_dist[1] ** 2)

        # for dynamic accelerations check if it's possible to fit 
        if accel_type == 'dynamic' or accel_type == 'destructive':
            old_accel = acceleration
            acceleration, deceleration = check_dynamic_acceleration(max(self.periods), dist, start_vel,
                                                                    max_vel, end_vel, acceleration,
                                                                    deceleration, acceleration/2)
            if acceleration == 0:
                self.da_revert += 1
                acceleration = old_accel
                accel_type = 'constant'
            if deceleration == 0:
                self.da_revert += 1
                deceleration = old_accel
                decel_type = 'constant'

        # vel, accel
        start_vel_xy = start_vel * x_ratio, start_vel * y_ratio
        x_vel = max_vel * x_ratio
        y_vel = max_vel * y_ratio
        x_end_vel = end_vel * x_ratio
        y_end_vel = end_vel * y_ratio
        x_accel = acceleration * x_ratio
        y_accel = acceleration * y_ratio
        x_decel = deceleration * x_ratio
        y_decel = deceleration * y_ratio

        # linear acceleration things
        accel_time = (max_vel-start_vel) / acceleration
        accel_dist = .5 * (max_vel + start_vel) * accel_time
        decel_time = (max_vel-end_vel) / deceleration
        decel_dist = .5 * (max_vel + end_vel) * decel_time

        # check validity of move
        if dist < accel_dist - self.eps:
            raise ValueError("can't fully accelerate. {} - {}".format(dist, accel_dist))
        # accel decel only move
        elif dist < (accel_dist + decel_dist) - self.eps:
            raise ValueError("can't fully accelerate and decelerate. {} - {} + {}".format(dist, accel_dist, decel_dist))

        # possible amount of cruising
        cruise_dist = dist - (accel_dist + decel_dist)
        if cruise_dist < 0:
            cruise_dist = 0
        cruise_time = cruise_dist / (max_vel)
        cruise_time_samples = int(round(cruise_time / (1/self.sample_rate)))

        data = {
            'accel_start': 0,
            'accel_end': 0,
            'decel_start': 0,
            'decel_end': 0
        }

        # accel->cruise->decel move
        iters = 0
        data['accel_start'] = 0
        # due to jerk, it's possible no accel is needed
        if accel_time > 0:
            move_time = self._accelerate_xy(start_vel_xy, (x_vel, y_vel), (x_accel, y_accel), accel_time, accel_type)
            iters += move_time
        data['accel_end'] = iters
        # cruise for some time period
        # make sure we're at right speed, can happen if no accel this move or decel last move
        cruise_vel = (x_vel, y_vel)
        for a_ind, a in enumerate(self.axes):
            a[0].velocity = (cruise_vel[a_ind], 0)
        for _ in range(cruise_time_samples):
            self._physics_update()
            iters += 1
        data['decel_start'] = iters
        # due to jerk, it's possible no decel is needed
        if decel_time > 0:
            move_time = self._accelerate_xy((x_vel, y_vel), (x_end_vel, y_end_vel), (-x_decel, -y_decel), decel_time, decel_type)
            iters += move_time
        data['decel_end'] = iters

        # there will always be some error due to simulation constraints check if any are very large
        # might be caused by bad acceleration function
        for a_ind, a in enumerate(self.axes):
            delta_pos = a[0].position[0] - start_pos[a_ind]
            error_pos = xy_dist[a_ind] - delta_pos
            if abs(error_pos) > self.eps:
                print("Move {}. Axis {} didn't move all the way. Error {:.4f}mm, check acceleration function.".format(xy_dist, a_ind, error_pos))

        return self.axes_metrics, data


    def _move_xy(self, xy_dist, acceleration, deceleration, accel_t, cruise_t, decel_t, start_vel, max_vel, end_vel, accel_type):
        # save start position for accuracy check at the end of this function
        start_pos = [a[0].position[0] for a in self.axes]
        # assumes max_vel is total equating to how much plastic can be extruded
        # decompose x/y amounts of move, accel, vels
        x_ratio, y_ratio = decompose_xy(xy_dist[0], xy_dist[1])
        accel_type, decel_type = accel_type, accel_type

        # dist
        dist = math.sqrt(xy_dist[0] ** 2 + xy_dist[1] ** 2)

        # vel, accel
        start_vel_xy = start_vel * x_ratio, start_vel * y_ratio
        x_vel = max_vel * x_ratio
        y_vel = max_vel * y_ratio
        x_end_vel = end_vel * x_ratio
        y_end_vel = end_vel * y_ratio
        x_accel = acceleration * x_ratio
        y_accel = acceleration * y_ratio
        x_decel = deceleration * x_ratio
        y_decel = deceleration * y_ratio

        cruise_time_samples = int(cruise_t * self.sample_rate)

        data = {
            'accel_start': 0,
            'accel_end': 0,
            'decel_start': 0,
            'decel_end': 0
        }

        # accel->cruise->decel move
        iters = 0
        data['accel_start'] = 0
        # due to jerk, it's possible no accel is needed
        if accel_t > 0:
            move_time = self._accelerate_xy(start_vel_xy, (x_vel, y_vel), (x_accel, y_accel), accel_t, accel_type)
            iters += move_time
        data['accel_end'] = iters
        # cruise for some time period
        # make sure we're at right speed, can happen if no accel this move or decel last move
        cruise_vel = (x_vel, y_vel)
        for a_ind, a in enumerate(self.axes):
            a[0].velocity = (cruise_vel[a_ind], 0)
        for _ in range(cruise_time_samples):
            self._physics_update()
            iters += 1
        data['decel_start'] = iters
        # due to jerk, it's possible no decel is needed
        if decel_t > 0:
            move_time = self._accelerate_xy((x_vel, y_vel), (x_end_vel, y_end_vel), (-x_decel, -y_decel), decel_t, decel_type)
            iters += move_time
        data['decel_end'] = iters

        # there will always be some error due to simulation constraints check if any are very large
        # might be caused by bad acceleration function
        for a_ind, a in enumerate(self.axes):
            delta_pos = a[0].position[0] - start_pos[a_ind]
            error_pos = xy_dist[a_ind] - delta_pos
            if abs(error_pos) > self.eps:
                print("Move {}. Axis {} didn't move all the way. Error {:.4f}mm, check acceleration function.".format(xy_dist, a_ind, error_pos))

        return self.axes_metrics, data

    def _accelerate_xy(self, start_vel, end_vel, accel, accel_time, accel_type):
        # assumes all params are (x, y) tuples except accel_time since it's the same
        iteration = 0
        # account for simulation sample rate
        accel = [a / self.sample_rate for a in accel]
        accel_samples = accel_time * self.sample_rate
        # handle instantaneous change in start vel if present
        for a_ind, a in enumerate(self.axes):
            a[0].velocity = (start_vel[a_ind], 0)

        while iteration / accel_samples < 1:
            for a_ind, a in enumerate(self.axes):
                driver = a[0]
                # uniform accel, dynamic accel is also constant
                if accel_type == 'constant' or accel_type == 'dynamic':
                    accel_this_step = accel[a_ind]
                # smoothstep velocity
                elif accel_type == 'smoothstep':
                    accel_this_step = smoothstep_accel(iteration/accel_samples) * accel[a_ind]
                # smootherstep velocity
                elif accel_type == 'smootherstep':
                    accel_this_step = smootherstep_accel(iteration/accel_samples) * accel[a_ind]
                elif accel_type == 'destructive':
                    accel_this_step = destruct(iteration, accel_samples) * accel[a_ind]
                else:
                    raise AttributeError('Invalid accel_type got {}'.format(accel_type))
                driver.velocity = (driver.velocity[0] + accel_this_step, 0)
            iteration += 1
            self._physics_update()

        # make sure any rounding errors are handled
        for a_ind, a in enumerate(self.axes):
            a[0].velocity = (end_vel[a_ind], 0)
        return iteration

    def _physics_update(self):
        # Update physics
        dt = 1/self.sample_rate
        self.space.step(dt)

        # save positions & accel of objects
        for ax, met in zip(self.axes, self.axes_metrics):
            driver, platform = ax
            met['expected'].append(driver.position[0] + self.rest_length)
            met['actual'].append(platform.position[0])
            met['error'].append(platform.position[0] - (driver.position[0] + self.rest_length))

    ### The rest of this is to create objects in PyMunk ###
    def create_poly(self, points, mass = 5.0, pos = (0,0)):
        moment = pm.moment_for_poly(mass, points)
        body = pm.Body(mass, moment, body_type=pm.Body.KINEMATIC)
        body.position = Vec2d(pos)
        shape = pm.Poly(body, points)
        shape.friction = 0.5
        shape.collision_type = self.COLLTYPE_DEFAULT
        self.space.add(body, shape)
        return body, shape

    def create_box(self, pos, size = 10, mass = 5.0):
        box_points = [(-size, -size), (-size, size), (size,size), (size, -size)]
        return self.create_poly(box_points, mass = mass, pos = pos)

    def create_circle(self, point, mass=1.0, radius=15.0):
        moment = pm.moment_for_circle(mass, 0.0, radius)
        ball_body = pm.Body(mass, moment)
        ball_body.position = Vec2d(point)

        ball_shape = pm.Circle(ball_body, radius)
        ball_shape.friction = 1.5
        ball_shape.collision_type = self.COLLTYPE_DEFAULT
        self.space.add(ball_body, ball_shape)
        return ball_body, ball_shape

    def create_axis(self, pos_y, stiffness, damping_coef):
        square, square_shape = self.create_box((-self.rest_length, pos_y))
        circle_mass = 1
        circle, circle_shape = self.create_circle((0, pos_y), mass=circle_mass, radius=10)

        # critical damping, best case
        critical_damping = 2*math.sqrt(stiffness*circle_mass)
        damping = critical_damping * damping_coef

        # spring between them
        spring = pm.constraint.DampedSpring(square, circle, [0, 0], [0, 0], self.rest_length, stiffness, damping)
        self.space.add(spring)
        return square, circle

