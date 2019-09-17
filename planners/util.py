import numpy as np

def euclidean_dist(xy_pair):
    x, y = xy_pair
    return np.sqrt(x**2 + y**2)

def stiffness_from_period_damping(period, damping):
    # assumes mass == 1
    # T from mass stiff -> T = 2*np.pi*np.sqrt(mass/stiffness)
    # solve for stiffness
    # stiffness_no_damping = 1 / (period / (2 * np.pi)) ** 2

    # however damping > 0 affects this into
    # we are only given damped period and estimate of damping
    # damped_period = (1-damping)**0.5 * natural_period
    # damped_period = (1-damping)**0.5 * (2*np.pi*np.sqrt(mass/stiffness))
    # 2*np.pi*np.sqrt(mass/stiffness) = damped_period / (1-damping)**0.5
    # np.sqrt(mass/stiffness) = (damped_period / (1-damping)**0.5) / (2*np.pi)
    # 1 / stiffness = ((damped_period / (1-damping)**0.5) / (2*np.pi))**2

    # just in case period is 0, it's never 0
    period = max(period, 0.0001)

    if damping == 1:
        return 1 / (period / (2 * np.pi)) ** 2

    return 1 / ((period / np.sqrt(1-damping)) / (2*np.pi))**2

def decompose_xy(x, y):
    """
    This code assumes x, y are offsets from current position (ie: distance on each axis)
    """
    d_total = euclidean_dist((x, y))
    # what do each axis contribute to the move
    if d_total > 0:
        x_ratio = x / d_total
        y_ratio = y / d_total
        return x_ratio, y_ratio
    return 1, 1

def angle_between(xy1, xy2):
    # from https://stackoverflow.com/a/13849249
    n1, n2 = np.linalg.norm(xy1), np.linalg.norm(xy2)
    u1 = [c / n1 for c in xy1]
    u2 = [c / n2 for c in xy2]
    return np.rad2deg(np.arccos(np.clip(np.dot(u1, u2), -1, 1)))

def jerk_calculate_exit_entry_vel(xy1, vel1, xy2, vel2, jerk):
    x_exit_vel, y_exit_vel = [d * vel1 for d in decompose_xy(*xy1)]
    x_entry_vel, y_entry_vel = [d * vel2 for d in decompose_xy(*xy2)]
    delta_x = abs(x_exit_vel-x_entry_vel)
    delta_y = abs(y_exit_vel-y_entry_vel)
    delta = max(delta_x, delta_y)
    # don't let speed increase
    scale = min(jerk / delta, 1)
    return scale

def jerk_calculate_junction(xy1, vel1, xy2, vel2, jerk):
    x_exit_vel, y_exit_vel = [d * vel1 for d in decompose_xy(*xy1)]
    x_entry_vel, y_entry_vel = [d * vel2 for d in decompose_xy(*xy2)]
    delta_x = abs(x_exit_vel-x_entry_vel)
    delta_y = abs(y_exit_vel-y_entry_vel)
    delta = max(delta_x, delta_y)
    return jerk / delta

def jerk_calculate(xy1, vel1, xy2, vel2):
    x_exit_vel, y_exit_vel = [d * vel1 for d in decompose_xy(*xy1)]
    x_entry_vel, y_entry_vel = [d * vel2 for d in decompose_xy(*xy2)]
    delta_x = abs(x_exit_vel-x_entry_vel)
    delta_y = abs(y_exit_vel-y_entry_vel)
    return max(delta_x, delta_y)

def find_max_vel(dist, start_vel, end_vel, acceleration):
    start_vel2, end_vel2 = start_vel ** 2, end_vel ** 2
    max_end_vel2 = start_vel**2 + 2*acceleration*dist
    max_start_vel2 = end_vel**2 + 2*acceleration*dist
    intersection_dist_vel = _line_intersection((0, start_vel2), (dist, max_end_vel2), (0, max_start_vel2), (dist, end_vel2))
    # move is impossible probably should alert user
    # if intersection_dist_vel[0] > dist or intersection_dist_vel[0] < 0:
        # print('this move is impossible', intersection_dist_vel)
    return np.sqrt(intersection_dist_vel[1])

def smoothstep_accel(x):
    return 6*x - 6*x**2

def smootherstep_accel(x):
    return 30*x**2 - 60*x**3 + 30*x**4

def smoothererstep_accel(x):
    return 140*x**3 - 420*x**4 + 420*x**5 - 140*x**6

def jerk_limited_accel_ramp(iteration, time):
    # TODO: this is a jerk limited accel ramp
    # each ramp is 1/4 of the total time
    # the sum of the ramps is 1/4*1/2
    # so the total sum we have is 0.75 = .125 + .5 + .125
    # we want accel_time to not change so scale it
    ramp_time = 0.25
    cruise_time = 1 - ramp_time*2
    scalar = 1/(cruise_time + ramp_time)

    accel_time = int(np.ceil(time))
    normalized_time = iteration / accel_time
    if normalized_time < ramp_time:
        scalar_time = normalized_time / ramp_time
        return scalar_time * scalar
    elif normalized_time < cruise_time + ramp_time:
        return scalar
    else:
        scalar_time = 1 - ((normalized_time - (cruise_time + ramp_time)) / ramp_time)
        return scalar_time * scalar

def check_dynamic_acceleration(period, dist, start_vel, max_vel, end_vel, acceleration, deceleration, min_accel):
    # linear acceleration things
    if max_vel - start_vel > 0:
        acceleration = max((max_vel-start_vel) / period, min_accel)
    accel_time = (max_vel-start_vel) / acceleration
    accel_dist = .5 * (max_vel + start_vel) * accel_time
    if max_vel - end_vel > 0:
        deceleration = max((max_vel-end_vel) / period, min_accel)
    decel_time = (max_vel-end_vel) / deceleration
    decel_dist = .5 * (max_vel + end_vel) * decel_time
    # not enough distance to accel
    if dist < accel_dist:
        return 0, 0
    # not enough distance to accel & decel
    elif dist < (accel_dist + decel_dist):
        return 0, 0
    return acceleration, deceleration

def destruct(iteration, accel_time):
    # dynamic signal canceling
    accel_time = accel_time
    pulse_length = accel_time * 0.416666
    offset = accel_time * 0.08888888888
    scale = 1.485714286
    # pulse_length = accel_time * 0.1166667
    # offset = accel_time * 0.1
    # scale = 2.2632653

    t1 = offset
    t2 = pulse_length
    t3 = accel_time - (pulse_length + offset)

    linear_dist = 0.5*accel_time**2

    # new signal is a piecewise linear fn
    v1 = offset
    v2 = v1+(scale*pulse_length)
    v3 = v2+(accel_time-(pulse_length+offset))
    new_dist = 0.5*(0+v1)*offset + 0.5*(v1+v2)*pulse_length + 0.5*(v2+v3)*(accel_time-(pulse_length+offset))
    dist_scale = linear_dist/new_dist
    if iteration < offset:
        return dist_scale
    elif iteration < offset + pulse_length:
        return scale*dist_scale
    else:
        return dist_scale

def _line_intersection(l1_p1, l1_p2, l2_p1, l2_p2):
    xdiff = (l1_p1[0] - l1_p2[0], l2_p1[0] - l2_p2[0])
    ydiff = (l1_p1[1] - l1_p2[1], l2_p1[1] - l2_p2[1]) #Typo was here

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(l1_p1, l1_p2), det(l2_p1, l2_p2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def gcode_read(gcode_file):
    with open(gcode_file, 'r') as in_file:
        gcode_moves = []
        lines = in_file.readlines()
        for l_ind, l in enumerate(lines):
            # not comment line
            if l[0] != ';':
                try:
                    split = l.split(' ')
                    gcommand = split[0]
                    # is a g1 command
                    if gcommand.lower() == 'g1' or gcommand.lower() == 'g0':
                        x, y = split[1], split[2]
                        # can be feed rate, if so x, y are 2, 3
                        if x[0].lower() == 'f':
                            if len(split) > 3:
                                x, y = split[2], split[3]
                            else:
                                continue
                        # make sure these are x & y
                        if x[0].lower() != 'x':
                            raise NotImplementedError('GCode Parse Error: line {}, does not have X'.format(l_ind))
                        elif y[0].lower() != 'y':
                            raise NotImplementedError('GCode Parse Error: line {}, does not have Y'.format(l_ind))
                        else:
                            gcode_moves.append((float(x[1:]), float(y[1:])))
                except Exception as e:
                    raise IOError("Unable to parse Gcode line {}. {}".format(l_ind, l))
    return gcode_moves


