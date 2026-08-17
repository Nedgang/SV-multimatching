#############
# FUNCTIONS #
#############
def list_merged_intervals(list_intervals: list((int, int))) -> list:
    """
    From a sorted list of intervals [(start, end), (start, end)], return a list where
    overlapping intervals are merged together.
    Input: [(int, int)]
    Output: [(int, int)]
    """
    if list_intervals == []:
        return []
    else:
        return_list = []
        interval_start = list_intervals[0][0]
        interval_end = list_intervals[0][1]
        for i in range(1, len(list_intervals)):
            next_start = list_intervals[i][0]
            next_end = list_intervals[i][1]
            if interval_end < next_start:
                return_list.append((interval_start, interval_end))
                interval_start = next_start
                interval_end = next_end
            elif interval_end < next_end:
                interval_end = next_end
            else:
                pass
        return_list.append((interval_start, interval_end))
    return return_list


def overlap_size(interval_a: (int, int), interval_b: (int, int)) -> int:
    """
    Return overlap size from 2 tuple intervals.
    Input: interval_a, interval_b, both (int, int)
    Output: int
    """
    if interval_a[1] < interval_b[0] or interval_a[0] > interval_b[1]:
        return 0
    else:
        # Need to add 1 to have exact overlap value (SNV case of start == end)
        return (
            min([interval_a[1], interval_b[1]])
            - max([interval_a[0], interval_b[0]])
            + 1
        )


def is_list_intervals_in_limits(
    list_intervals: list((int, int)), limit: int, var_start: int, var_end: int
) -> bool:
    """
    Check if there is something in the intervals list, and if the extremities are in
    the acceptable threshold.
    """
    if list_intervals == []:
        return False
    elif (
        list_intervals[0][0] > var_start - limit
        and list_intervals[0][0] < var_start + limit
        and list_intervals[-1][1] > var_end - limit
        and list_intervals[-1][1] < var_end + limit
    ):
        return True
    else:
        return False
