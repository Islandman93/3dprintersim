import numpy as np
from .util import decompose_xy, jerk_calculate_junction, find_max_vel, euclidean_dist, jerk_calculate_exit_entry_vel


class MarlinPlanner:
    def __init__(self, acceleration, max_vel, jerk):
        self.acceleration = acceleration
        # assumes max_vel is total equating to how much plastic can be extruded
        self.max_vel = max_vel
        self.jerk = jerk
        self.moves = []
        self._last_velocity = 0

    def flush(self):
        reachable_vels = []
        # iterate backwards assuming last move end_vel = 0.05
        end_vel = 0.05
        for i in range(len(self.moves)-1, -1, -1):
            move = self.moves[i]
            dist = np.sqrt(move[0] ** 2 + move[1] ** 2)
            # maximum starting velocity is full decel to current end_vel
            max_start_vel = np.sqrt(end_vel ** 2 + (2*self.acceleration*dist))
            # TODO: this has an actual function I just don't know what right now
            if i != 0:
                jerk_allowed = self.jerk * jerk_calculate_junction(self.moves[i-1], self.jerk, move, self.jerk, self.jerk)
            else:
                jerk_allowed = self.jerk / 2
            # if max_start_vel > jerk allowed use jerk as start_vel
            if max_start_vel >= jerk_allowed:
                start_vel = jerk_allowed
                max_end_vel = np.sqrt(start_vel ** 2 + (2*self.acceleration*dist))
                # requested start vel to end vel may be impossible
                # TODO: refactor this into function
                if end_vel > max_end_vel:
                    # iterate forward from this point based on max end vel until it can be reached
                    for j in range(len(reachable_vels)-1, -1, -1):
                        next_move = reachable_vels[j]
                        next_ind = i + j + 1
                        next_dist = euclidean_dist(self.moves[next_ind])
                        # can't reach requested env vel update and iterate
                        if max_end_vel < next_move[0]:
                            new_reachable_end_vel = np.sqrt(max_end_vel ** 2 + (2*self.acceleration*next_dist))
                            # done just adjust cruise vel
                            if new_reachable_end_vel > reachable_vels[j][-1]:
                                cruise_vel = min(find_max_vel(next_dist, max_end_vel, reachable_vels[j][-1], self.acceleration), self.max_vel)
                                reachable_vels[j] = (max_end_vel, cruise_vel, reachable_vels[j][-1])
                                break
                            # otherwise need to iterate
                            else:
                                cruise_vel = min(find_max_vel(next_dist, max_end_vel, new_reachable_end_vel, self.acceleration), self.max_vel)
                                reachable_vels[j] = (max_end_vel, cruise_vel, new_reachable_end_vel)
                                max_end_vel = new_reachable_end_vel
                        # end_vel is reachable
                        else:
                            break
                    # finally set end vel for this move
                    end_vel = np.sqrt(start_vel ** 2 + (2*self.acceleration*dist))
                    cruise_vel = end_vel
                else:
                    cruise_vel = min(find_max_vel(dist, start_vel, end_vel, self.acceleration), self.max_vel)
            else:
                start_vel = min(max_start_vel, self.max_vel)
                cruise_vel = start_vel
            reachable_vels.append((start_vel, cruise_vel, end_vel))
            end_vel = start_vel

        # iterate forwards assuming first move start_vel = jerk / 2
        # take the minimum of the reachable velocities of both passes
        vels = reachable_vels[::-1]
        return vels

    def add(self, xy_dist):
        self.moves.append(xy_dist)

    def add_move(self, xy_dist, next_xy_dist=None):
        # import pudb; pudb.set_trace()
        # decompose x/y amounts of move, accel, vels
        x_ratio, y_ratio = decompose_xy(xy_dist[0], xy_dist[1])

        # total dist
        dist = np.sqrt(xy_dist[0] ** 2 + xy_dist[1] ** 2)
        if next_xy_dist is not None:
            x_ratio_next, y_ratio_next = decompose_xy(next_xy_dist[0], next_xy_dist[1])
            dist_next = np.sqrt(next_xy_dist[0] ** 2 + next_xy_dist[1] ** 2)
        else:
            dist_next = 0

        # starting velocity is last move's end vel
        if self._last_velocity == 0:
            # marlin planner uses half jerk for moves starting at 0
            start_vel = self.jerk / 2
        else:
            start_vel = self._last_velocity

        # estimate the lowest cruise velocity of the this move and the next
        max_vel = find_max_vel(dist, start_vel, 0, self.acceleration)
        max_vel = min(max_vel, self.max_vel)
        if dist_next != 0:
            max_vel_next = find_max_vel(dist_next, 0, 0, self.acceleration)
            max_vel_next = min(max_vel_next, self.max_vel)
        else:
            max_vel_next = 0

        if dist_next != 0:
            # check what we need to decel to 
            vel_scale = jerk_calculate_exit_entry_vel(xy_dist, max_vel, next_xy_dist, max_vel_next, self.jerk)
            end_vel = max_vel * vel_scale
        else:
            # marlin decelerate to 0.05
            end_vel = 0.05
        self._last_velocity = end_vel

        return start_vel, max_vel, end_vel

    def get_next_move(self):
        pass


class Move:
    def __init__(self, start_vel, max_vel, end_vel, acceleration):
        pass
