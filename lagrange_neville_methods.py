
def validate_and_prepare_points(points_table, x_to_find):
    sorted_points = sorted(points_table)

    for point in sorted_points:
        if point[0] == x_to_find:
            return None, None, point[1]

    if x_to_find < sorted_points[0][0] or x_to_find > sorted_points[-1][0]:
        return None, None, "This is extrapolation, not interpolation"

    x_list = []
    y_list = []

    for point in sorted_points:
        x_list.append(point[0])
        y_list.append(point[1])

    return x_list, y_list, None

def lagrange_method(points_table, x_to_find):
    x_list, y_list, result = validate_and_prepare_points(points_table, x_to_find)
    if result is not None:
        return result

    result = 0
    for i in range(len(x_list)):
        result += y_list[i] * Li_x(x_list, i, x_to_find)

    return result

def Li_x(x_list,i,x):
    res1=1
    res2=1
    for j in range(len(x_list)):
        if i==j:
            continue
        res1 *=(x-x_list[j])
        res2 *=(x_list[i]-x_list[j])
    num = res1/res2
    return num

def neville_method(points_table, x_to_find):
    x_list, y_list, result = validate_and_prepare_points(points_table, x_to_find)
    if result is not None:
        return result

    memo = {}

    def neville_recursive(i, j):
        if (i, j) in memo:
            return memo[(i, j)]

        if i == j:
            memo[(i, j)] = y_list[i]
            return memo[(i, j)]

        left_value = neville_recursive(i, j - 1)
        right_value = neville_recursive(i + 1, j)

        result = ((x_to_find - x_list[j]) * left_value +
                  (x_list[i] - x_to_find) * right_value) / (x_list[i] - x_list[j])

        memo[(i, j)] = result
        return result

    return neville_recursive(0, len(x_list) - 1)


def main():
    points_table = [
        [1, 0],
        [1.2, 0.112463],
        [1.3, 0.167996],
        [1.4, 0.222709]
    ]

    x_to_find = 1.28

    lagrange_result = lagrange_method(points_table, x_to_find)
    neville_result = neville_method(points_table, x_to_find)

    print("Points table:", points_table)
    print("X value:", x_to_find)
    print("lagrange result:", lagrange_result)
    print("neville_method:", neville_result)



if __name__ == "__main__":
    main()
