#! /usr/bin/python3.3
import copy


def Zeidel_iteration(matrix, rightCol, inputVector):
    result = [None] * len(inputVector)
    for i,k in enumerate(rightCol):
        result[i] = rightCol[i]
        for j,k in enumerate(rightCol):
            result[i] += matrix[i][j] * inputVector[j] if i <= j else matrix[i][j] * result[j] 
    return result


def Zeidel_met(matrix, rightCol, vector, epselon):    
    readyMatrix, readyCol = prepare_matrix(matrix, rightCol)
    i=0
    while True:
        vectorNew = Zeidel_iteration(readyMatrix, readyCol, vector)
        i += 1
        if exit_test(vector,vectorNew) < epselon:
            break
        vector = vectorNew
    return [round(x,4) for x in vector]


def exit_test(vector1,vector2):
    return max([abs(k-v) for k,v in zip(vector1,vector2)])
           

def prepare_matrix(matrix, rightCol):
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i != j:
                matrix[i][j] /= -matrix[i][i]
        rightCol[i] /= matrix[i][i]
        matrix[i][i] = 0
    return matrix,rightCol
