from gauss_seidel import gauss_seidel


def linear_interpolation(points, x_to_find):

    sorted_points = sorted(points, key=lambda point: point[0])

    for i in range(len(sorted_points) - 1):
        x1, y1 = sorted_points[i]
        x2, y2 = sorted_points[i + 1]

        if x1 <= x_to_find <= x2:
            if x1 == x2:
                raise ValueError("Cannot interpolate between two points with the same x value.")

            return y1 + (x_to_find - x1) * (y2 - y1) / (x2 - x1)

    raise ValueError("The x value is outside the range of the points table.")


def polynomial_interpolation(points, x):
    points = sorted(points, key=lambda point: point[0])
    x_list = [point[0] for point in points]
    y_list = [point[1] for point in points]
    n = len(points)

    if n < 2:
        print("The points table must contain at least two points.")
        return

    for i in range(n):
        if x == x_list[i]:
            return y_list[i]

    is_between_points = False
    for i in range(n - 1):
        if x_list[i] < x < x_list[i + 1]:
            is_between_points = True
            break

    if not is_between_points:
        print("The x value is not between two x values in the table.")
        return

    matrix = []

    for i in range(n):
        row = []
        for j in range(n):
            row.append(x_list[i] ** j)
        matrix.append(row)

    results = []
    for y in y_list:
        results.append([y])

    b = gauss_seidel(matrix, results)

    p_x = 0

    for i in range(n):
        p_x = p_x + b[i] * (x ** i)

    return p_x


def main():
    points_table = [
        [0, 0],
        [1, 0.8415],
        [2, 0.9093],
        [3, 0.1411],
        [4, -0.7568],
        [5, -0.9589],
        [6, -0.2794]
    ]

    x_to_find = 2.5

    linear_result = linear_interpolation(points_table, x_to_find)
    polynomial_result = polynomial_interpolation(points_table, x_to_find)

    print("Points table:", points_table)
    print("X value:", x_to_find)
    print("Linear interpolation result:", linear_result)
    print("Polynomial interpolation result:", polynomial_result)


if __name__ == "__main__":
    main()
