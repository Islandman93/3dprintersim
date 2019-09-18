import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
import numpy as np

from simulation import Simulation
from planners.util import gcode_read, jerk_calculate, angle_between


EPSILON = 0.0125
parser = argparse.ArgumentParser()
subparser = parser.add_subparsers(dest='planner')

# default args for both planners
def add_default_args(parser):
    parser.add_argument('gcode', type=str, help='Location of GCode file, recommended to only have 1 layer')
    parser.add_argument('velocity', type=int, help='Maximum velocity')
    parser.add_argument('acceleration', type=int, help='Acceleration')
    parser.add_argument('period_x', type=float, help='Determines stiffness of the x spring in simulation')
    parser.add_argument('period_y', type=float, help='Determines stiffness of the y spring in simulation')
    parser.add_argument('-at', '--acceleration_type', type=str.lower, choices=
                        ['constant', 'smoothstep', 'smootherstep', 'dynamic', 'destructive',
                         'jerk-limited-smoothstep', 'jerk-limited-smootherstep'],
                        help='Acceleration type', default='constant')
    parser.add_argument('-o', '--origin', nargs='+', type=float, help='XY origin point', default=(0, 0))
    parser.add_argument('-sr', '--sample_rate', type=int, help='The number of physics steps in a second', default=40000)

# klipper planner args
parser_klipper = subparser.add_parser('klipper', help='Klipper planner')
add_default_args(parser_klipper)
parser_klipper.add_argument('square_corner_velocity', type=float, help='Square Corner Velocity')
parser_klipper.add_argument('-a2d', '--accel_to_decel', type=int, help='Ghost deceleration to provide "smoothed look-ahead"')

# marlin planner args
parser_marlin = subparser.add_parser('marlin', help='Marlin-style "jerk" limited planner')
add_default_args(parser_marlin)
parser_marlin.add_argument('jerk', type=float, help='Jerk as defined in marlin, the maximum instantaneous change in velocity')
parser_marlin.add_argument('--naive', action='store_true', help='Uses a naive planner that performs no move checking. Will not work in all cases')

args = parser.parse_args()

# parse gcode to moves
gcode_moves = gcode_read(args.gcode)
# convert to x, y offsets this is how the simulation processes moves
gcode_move_offsets = []
tiny_loc = []
num_tiny = 0
# find small moves
last = args.origin
first_move_position = args.origin
for m in gcode_moves:
    delta = [a - b for a, b in zip(m, last)]
    if np.sqrt(delta[0]**2 + delta[1]**2) < 0.0125 * 5:
        angle = angle_between(last_delta, delta)
        if angle < 1:
            print('{} Tiny move in a straight line from {} to {}. Angle {:.3f}'.format(len(gcode_move_offsets), last, m, angle_between(last_delta, delta)))
        tiny_loc.append(last)
        tiny_loc.append(m)
        num_tiny += 1

    gcode_move_offsets.append(delta)
    last = m
    last_delta = delta

if num_tiny > 0:
    print('GCode contains {} tiny moves < {:.3f}mm'.format(num_tiny, 0.0125 * 5))

# init Simulation
sim = Simulation(2, (args.period_x, args.period_y), sample_rate=args.sample_rate)

# klipper planner with square_corner_velocity
if args.planner == 'klipper':
    from planners.klipperplanner import KlipperPlanner

    # only used for jerk-limited acceleration profiles
    jerk = 0.6 * args.acceleration * (1 / max(args.period_x, args.period_y))
    klip_plan = KlipperPlanner(args.acceleration, args.velocity, args.square_corner_velocity, args.accel_to_decel, args.acceleration_type, jerk=jerk)

    # add to planner's move queue
    current_pos = first_move_position
    for m_ind, m in enumerate(gcode_moves):
        try:
            klip_plan.add_move(current_pos, m)
            current_pos = m
        except Exception as e:
            print(gcode_moves[m_ind-1])
            print('Error on move {}: {}'.format(m_ind, m))
            print(gcode_moves[m_ind+1])
            raise e

    # get planned moves start, cruise, end velocities
    planned_moves = klip_plan.flush()

    # klipper returns a list of move objects
    # simulation needs move offsets
    for move_ind, (move_offset, planned_move) in enumerate(zip(gcode_move_offsets, planned_moves)):
        try:
            sim._move_xy(move_offset, planned_move.accel, planned_move.decel, planned_move.accel_t, planned_move.cruise_t, planned_move.decel_t, planned_move.start_v, planned_move.cruise_v,
                        planned_move.end_v, args.acceleration_type)
        except Exception as e:
            print('Error on move {}'.format(move_ind))
            raise e

# marlin planner "jerk" limited
else:
    from planners.marlinplanner import MarlinPlanner
    mar_plan = MarlinPlanner(args.acceleration, args.velocity, args.jerk)

    # marlin planner works in offset distances
    for i in range(len(gcode_move_offsets)):
        cur_move = gcode_move_offsets[i]
        if i == len(gcode_move_offsets) - 1:
            next_dist = None
        else:
            next_dist = gcode_move_offsets[i+1]

        # old code that does a 1 move lookahead only based on lowest bound of acheivable velocity
        if args.naive:
            sv, mv, ev = mar_plan.add_move(cur_move, next_dist)
            sim.move_xy(cur_move, args.acceleration, sv, mv, ev, args.acceleration_type)
        else:
            # add to planner
            mar_plan.add(cur_move)

    # if not naive planner, get the adjusted velocities from planner
    if not args.naive:
        planned_moves = mar_plan.flush()
        last_move = (0, 0), (0, 0, 0)
        for move_ind, (move_offset, planned_vels) in enumerate(zip(gcode_move_offsets, planned_moves)):
            # temp code to make sure the planner isn't exceeding jerk
            real_jerk = jerk_calculate(last_move[0], last_move[1][-1], move_offset, planned_vels[0])
            real_jerk = 0 if np.isnan(real_jerk) else real_jerk
            if real_jerk > args.jerk + EPSILON:
                print('move {} exceeds jerk'.format(move_ind), real_jerk)

            # send move to simulator
            sim.move_xy(move_offset, args.acceleration, planned_vels[0], planned_vels[1], planned_vels[2], args.acceleration_type)
            last_move = move_offset, planned_vels

# wait for springs to settle for .05 second
# for i in range(int(0.05*args.sample_rate)):
    # sim._physics_update()

SAMPLE_RATE = args.sample_rate
# 2 axis metrics
x, y = sim.axes_metrics
total_len = len(x['actual'])
total_error = np.sum(np.abs(x['error'])) + np.sum(np.abs(y['error']))
total_time = total_len / SAMPLE_RATE
print('Error: {:.4f}mm/s | Time: {:.3f}s'.format(total_error/SAMPLE_RATE, total_time))
if args.acceleration_type == 'dynamic' or args.acceleration_type == 'destructive':
    print('Dynamic accel number of revertions. {}/{} = {:.2f}%'.format(sim.da_revert, len(gcode_moves)*2, sim.da_revert/(len(gcode_moves)*2)))
# plot results with subsampling so it will actually render
# generally 1000 points per second is enough 
num_samples = 1000 * total_time
samples = np.sort(np.random.choice(total_len, int(num_samples), replace=False))

def sample(array):
    return np.asarray(array)[samples]

fig = plt.figure(constrained_layout=True)
gspec = gridspec.GridSpec(2, 2, fig)
# one full horizontal plot
ax = fig.add_subplot(gspec[:, 0])
ax.set_title('position')
# the simulation isn't perfect due to rounding errors, plot the exact gcode in black
# this should also visually check planner errors
# in case it gets really off
# subtract first move to put at origin
gcode_array = np.asarray(gcode_moves)
tiny_loc = np.asarray(tiny_loc)
gcode_array[:, 0] -= first_move_position[0]
gcode_array[:, 1] -= first_move_position[1]
gcode_line = ax.plot(gcode_array[:, 0], gcode_array[:, 1], 'k', label='gcode')
# if tiny moves
if num_tiny > 0:
    tiny_loc[:, 0] -= first_move_position[0]
    tiny_loc[:, 1] -= first_move_position[1]
    ax.plot(tiny_loc[:, 0], tiny_loc[:, 1], 'rx', label='tiny_move')
# compare driver and actual position
expected_pos_line, = ax.plot(sample(x['expected']), sample(y['expected']), 'orange', label='expected')
actual_pos_line, = ax.plot(sample(x['actual']), sample(y['actual']), label='actual')
ax.legend()
# two lower plots
ax = fig.add_subplot(gspec[0, 1])
ax.set_title('error')
x_error_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(x['error']), label='x')
y_error_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(y['error']), label='y')
ax.legend()
# plt.plot(np.sqrt(np.asarray(y['error'])**2+np.asarray(x['error'])**2), label='total')
ax = fig.add_subplot(gspec[1, 1])
ax.set_title('vel')
x_vel_expec, y_vel_expec = np.diff(x['expected']) * SAMPLE_RATE, np.diff(y['expected']) * SAMPLE_RATE
x_vel_actual, y_vel_actual = np.diff(x['actual']) * SAMPLE_RATE, np.diff(y['actual']) * SAMPLE_RATE
# toss on a zero since we are diffing
x_vel_actual = np.concatenate(([0], x_vel_actual))
y_vel_actual = np.concatenate(([0], y_vel_actual))
x_vel_expec = np.concatenate(([0], x_vel_expec))
y_vel_expec = np.concatenate(([0], y_vel_expec))
x_vel_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(x_vel_actual), label='x')
y_vel_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(y_vel_actual), label='y')
x_vel_expec_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(x_vel_expec), '-', color='teal', label='x_expected')
y_vel_expec_line, = ax.plot(sample(np.arange(total_len) / SAMPLE_RATE), sample(y_vel_expec), '-', color='orange', label='y_expected')
ax.legend()
plt.show()

# slice a certain number to show at each step
# show_num = 100
# slices = []
# for i in range(int(len(x['expected'])//show_num)):
    # slices.append(i * show_num)
# slices.append(len(x['expected']))

# def update(num):
    # time = np.arange(num) / SAMPLE_RATE
    # expected_pos_line.set_data(x['expected'][:num], y['expected'][:num])
    # actual_pos_line.set_data(x['actual'][:num], y['actual'][:num])
    # y_error_line.set_data(time, y['error'][:num])
    # x_error_line.set_data(time, x['error'][:num])
    # x_vel_line.set_data(time, x_vel_actual[:num])
    # y_vel_line.set_data(time, y_vel_actual[:num])


# anim = FuncAnimation(fig, update, interval=1, frames=slices, repeat=False)
# anim.save('simanim.gif', writer='imagemagick', fps=20)

