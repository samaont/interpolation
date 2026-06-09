def diagonal_dominate_check(matrix):
    n=len(matrix)
    for i in range(n):
        diagonal = abs(matrix[i][i])
        row_sum = sum(abs(matrix[i][j]) for j in range(n) if j!=i )
        if diagonal <= row_sum:
            return False
    return True

def make_diagonal_dominant(matrix,vector):
    n = len(matrix)
    matrix =[row[:] for row in matrix]
    vector =[row[:] for row in vector]
    for col in range(n):
        max_row = col
        for row in range(col+1,n):
            if abs(matrix[row][col]) > abs(matrix[max_row][col]):
                max_row = row
        matrix[col],matrix[max_row] = matrix[max_row],matrix[col]
        vector[col],vector[max_row] = vector[max_row],vector[col]
    if diagonal_dominate_check(matrix):
        return matrix,vector
    return None,None

def matrix_prepare(matrix,vector):
    if diagonal_dominate_check(matrix):
        return matrix,vector,True
    new_matrix,new_vector = make_diagonal_dominant(matrix,vector)
    if new_matrix is not None:
        print("matrix achieve diagonal dominate")
        return new_matrix,new_vector,True
    print("cant achieve diagonal dominate")
    return matrix,vector,False

def jacobi(matrixA,vectorB,tol=0.0001,max_iter=100):
    matrixA ,vectorB,flag= matrix_prepare(matrixA,vectorB)
    n=len(matrixA)
    x_old= [0.0]*n
    for iter in range(max_iter):
        x_new = [0.0]*n
        for i in range (n):
            sigma = sum(matrixA[i][j]*x_old[j] for j in range(n) if j!=i)
            x_new[i] = (vectorB[i][0]-sigma)/matrixA[i][i]

        if max(abs(x_new[i]-x_old[i]) for i in range(n)) < tol:
            if flag == False:
                print(f"Although there is no dominant diagonal converged after {iter} iteration ")
            else:
                print(f"converged after {iter} iteration ")
            return x_new

        x_old = x_new
    print("not converged")
    return x_old


def gauss_seidel(matrixA,vectorB,tol=0.0001,max_iter=100):
    matrixA,vectorB,flag = matrix_prepare(matrixA,vectorB)
    n = len(matrixA)
    x = [0.0] * n

    for iter in range(max_iter):
        x_old = x[:]
        for i in range(n):
            sigma = sum(matrixA[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (vectorB[i][0] - sigma) / matrixA[i][i]

        if max(abs(x[i] - x_old[i]) for i in range(n)) < tol:
            if flag == False:
                print(f"Although there is no dominant diagonal converged after {iter} iteration ")
            else:
                print(f"converged after {iter} iteration ")
            return x

    print("not converged")
    return x

def main():
    matrixA = [[4,2,0],[2,10,4],[0,4,5]]
    vectorB=[[2],[6],[5]]
    print("choose method:")
    print("1- jacobi")
    print("2- gauss seidel")
    choice=input("Enter your choice:")
    if choice=="1":
        jacobi(matrixA,vectorB)
    elif choice=="2":
        gauss_seidel(matrixA,vectorB)
    else:
        print("invalid choice")


if __name__ == "__main__":
    main()