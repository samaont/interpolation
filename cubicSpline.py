from gauss_seidel import gauss_seidel

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


def calculate_d_value(x_list, y_list, h_list, i):
    h_prev = h_list[i - 1]
    h_next = h_list[i]
    h_sum = h_prev + h_next

    return 6 / h_sum * (
        (y_list[i + 1] - y_list[i]) / h_next
        -
        (y_list[i] - y_list[i - 1]) / h_prev
    )


def build_internal_matrix_row(row_size, i, h_list, first_m_index):
    h_prev = h_list[i - 1]
    h_next = h_list[i]
    h_sum = h_prev + h_next

    row = [0.0] * row_size
    values = [
        (i - 1, h_prev / h_sum),
        (i, 2.0),
        (i + 1, h_next / h_sum)
    ]

    for m_index, value in values:
        row_index = m_index - first_m_index
        if 0 <= row_index < row_size:
            row[row_index] = value

    return row


def find_section_index(x_list, x_to_find):
    for i in range(len(x_list) - 1):
        if x_list[i] <= x_to_find <= x_list[i + 1]:
            return i
    return None


def calculate_spline_value(x_list, y_list, h_list, m_list, x_to_find, section_index):
    h = h_list[section_index]
    x_i = x_list[section_index]
    x_next = x_list[section_index + 1]
    y_i = y_list[section_index]
    y_next = y_list[section_index + 1]
    m_i = m_list[section_index]
    m_next = m_list[section_index + 1]

    return (
        m_i * (x_next - x_to_find) ** 3 / (6 * h)
        + m_next * (x_to_find - x_i) ** 3 / (6 * h)
        + (y_i - m_i * h ** 2 / 6) * (x_next - x_to_find) / h
        + (y_next - m_next * h ** 2 / 6) * (x_to_find - x_i) / h
    )


def cubic_spline_method(points_table, x_to_find, fTag0=None, fTagN=None):
    x_list, y_list, result = validate_and_prepare_points(points_table, x_to_find)
    if result is not None:
        return result

    h_list = []
    for i in range(len(x_list) - 1):
        h_list.append(x_list[i + 1] - x_list[i])

    internal_count = len(x_list) - 2

    if internal_count == 0:
        h = h_list[0]
        return y_list[0] + (y_list[1] - y_list[0]) * (x_to_find - x_list[0]) / h

    matrixA = []
    vectorB = []

    for i in range(1, len(x_list) - 1):
        matrixA.append(build_internal_matrix_row(internal_count, i, h_list, 1))
        vectorB.append([calculate_d_value(x_list, y_list, h_list, i)])

    internal_m_list = gauss_seidel(matrixA, vectorB)
    m_list = [0.0] + internal_m_list + [0.0]

    section_index = find_section_index(x_list, x_to_find)
    spline_value = calculate_spline_value(x_list, y_list, h_list, m_list, x_to_find, section_index)

    if fTag0 is None or fTagN is None:
        return spline_value

    full_matrixA = []
    full_vectorB = []

    first_row = [0.0] * len(x_list)
    first_row[0] = 2.0
    first_row[1] = 1.0
    first_d = 6 / h_list[0] * ((y_list[1] - y_list[0]) / h_list[0] - fTag0)
    full_matrixA.append(first_row)
    full_vectorB.append([first_d])

    for i in range(1, len(x_list) - 1):
        full_matrixA.append(build_internal_matrix_row(len(x_list), i, h_list, 0))
        full_vectorB.append([calculate_d_value(x_list, y_list, h_list, i)])

    last_row = [0.0] * len(x_list)
    last_row[-2] = 1.0
    last_row[-1] = 2.0
    last_d = 6 / h_list[-1] * (fTagN - (y_list[-1] - y_list[-2]) / h_list[-1])
    full_matrixA.append(last_row)
    full_vectorB.append([last_d])

    full_m_list = gauss_seidel(full_matrixA, full_vectorB)
    full_spline_value = calculate_spline_value(x_list, y_list, h_list, full_m_list, x_to_find, section_index)

    return spline_value, full_spline_value

def main():
    points_table = [
        [1, 0],
        [1.3, 0.6200],
        [1.6, 0.4554],
        [1.9, 0.2818],
        [2.2,0.1103]
    ]

    x_to_find = 1.5
    fTag0 = 1.0
    fTagN = -1.0

    result = cubic_spline_method(points_table, x_to_find, fTag0, fTagN)

    print(result)




if __name__ == "__main__":
    main()
