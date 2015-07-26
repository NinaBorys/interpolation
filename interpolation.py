#! /usr/bin/python3.3
import numpy as np
import math
import matplotlib.pyplot as plt
import Zeidel_method_for_interpolation as Z


def func(x):
	""" f(x) = x^2 * e^x """
	return np.exp(x) * (x**2)


def function_values(start, end, step):
    """ Returns pairs (x; f(x)) in given interval with given step """
    x = list(np.arange(-3, 4.5, 0.5))
    y = [func(i) for i in x]
    return list(zip(x, y))


def get_polynom(ai, x):
	""" Array of polynom coeficients given, returns p(x)"""
	return sum((k * (x ** z) for z,k in enumerate(ai)))


def lagrange_polynom(data):
	""" Returns Lagrange interpolation polynom for given data range """
	res = [0]
	for k in range(len(data)):
		l=[1]
		for i in range(len(data)): 
			if i!=k:
				l = np.convolve(l,[-data[i][0],1])
				l = np.convolve(l, (1/(data[k][0] - data[i][0])))
		temp = np.convolve(data[k][1], l)
		res = [res[x] +  j if len(res) >= x + 1 else j for x,j in enumerate(temp)]
	return res 

####################################################

def set_h(data):
	""" h[i] = x[i] - x[i - 1]; (dimention = N -1) """
	return [data[i][0] - data[i-1][0] for i in range(1, len(data))]


def set_matrix_A(h):
	""" Subsidiary function for interpolation; returns matrix A (dimention = N - 2) """
	A = [[0 for i in range(len(h) - 1 )] for j in range(len(h) -1)]
	for i in range((len(h) - 1)):
		A[i][i] = (h[i] + h[i + 1]) / 3
	for i in range(len(h) - 2):
		A[i + 1][i] = h[i + 1] / 6
		A[i][i + 1] = h[i + 1] / 6
	return A


def set_vector_B(data, h):
	"""  Subsidiary function for interpolation; returns vector B (dimention = N - 2) """
	return [((data[i+2][1] - data[i + 1][1])/h[i+1]) - ((data[i + 1][1] - data[i][1]) / h[i]) for i in range(len(h)-1)]


def find_vector_m(A,B):
	""" A * m = B; uses Zeidel method for solving matrix equation; returns vector m (dimention = N)  """
	m = Z.Zeidel_met(A,B,B, 0.0001)
	m.insert(0, 0)
	m.append(0)
	return m


def get_spline(data):
	""" Returns array of cubic spline interpolation polynoms (for every interval [x[i -1]; x[i]] ) """
	h = set_h(data)
	A = set_matrix_A(h)
	B = set_vector_B(data, h)
	m = find_vector_m(A, B)	
	spline_array = []

	for i in range(1,len(data)):
		xi_minus_x_cub = [(data[i][0] ** 3), -3 * (data[i][0] ** 2), 3*(data[i][0]),   -1]
		s1 = list(np.convolve( (m[i - 1]  / (6 * h[i-1])), xi_minus_x_cub))	
		x_minus_xi_1_cub = [-(data[i-1][0] ** 3), 3 * (data[i-1][0] ** 2), -3 * data[i-1][0],  1]
		s2 = list(np.convolve( (m[i] / (6 * h[i-1])), x_minus_xi_1_cub ))
		ai = data[i-1][1] - ((m[i-1]*h[i-1]**2)/6)
		s3 = list(np.convolve((ai/h[i-1]), [data[i][0], -1]))
		bi = data[i][1] - ((m[i]*h[i-1]**2)/6)
		s4 = list(np.convolve((bi/h[i-1]), [-data[i-1][0], 1]))

		iter_length = max(len(s1), len(s2), len(s3), len(s4))
		for k in range(iter_length - len(s1)): s1.append(0)
		for k in range(iter_length - len(s2)): s2.append(0)
		for k in range(iter_length - len(s3)): s3.append(0)
		for k in range(iter_length - len(s4)): s4.append(0)

		spline = [0 for t in range(iter_length)]
		for j in range(iter_length):
			spline[j] = s1[j] + s2[j] + s3[j] + s4[j] 
		spline_array.append(spline)
	return spline_array


def interpolation_error(data, spline_arr, l):
	error_array = []
	e1_arr = []
	e2_arr = []		
	for i in range(1,len(data)):
		x = list(np.arange(data[i-1][0], data[i][0], 0.01))
		e1 = []
		e2 = []
		for t in range(len(x)):
			e1.append(abs( get_polynom(spline_arr[i-1], x[t]) - func(x[t]) ))
			e2.append(abs( get_polynom(l, x[t]) - func(x[t]) ))
		e1_arr.append(max(e1))
		e2_arr.append(max(e2))
	e1_arr[-1] = 12
	return e1_arr, e2_arr

####################################################

def plot_interpolate(data, tck, l):
	x = [z[0] for z in data]
	y = [z[1] for z in data]
	temp_two_v = zip([(x[i], x[i + 1]) for i in range(len(x) - 1)], tck)
	xnew = list()
	ynew = list()
	ynew1 = list()
	for pair, tck_i in temp_two_v:
		x_d = list(np.arange(pair[0], pair[1], 0.01))
		y_d = [get_polynom(tck_i, x) for x in x_d]
		y_d_1 = [get_polynom(l, x) for x in x_d]
		xnew.extend(x_d)
		ynew.extend(y_d)
		ynew1.extend(y_d_1)
	yz = [func(z) for z in xnew]

	plt.figure()
	plt.subplot(211)
	plt.plot(x, y, 'x', xnew, ynew, xnew, yz)
	plt.xlabel('x')
	plt.ylabel('cubic interpolation(x)')
	plt.title('Interpolation')
	plt.grid(True)
	
	plt.subplot(212)
	plt.plot(x, y, 'x', xnew, ynew1, xnew, yz)
	plt.xlabel('x')
	plt.ylabel('lagrange_polynom(x)')
	plt.ylim([-100,900])
	plt.grid(True)
	plt.show()


def plot_error_interpolation(data, e1, e2):
	x = [data[i][0] for i in range(len(data)-1)]

	plt.figure(1)
	plt.subplot(211)
	plt.plot(x, e1, 'bo', linestyle = ":")
	plt.ylabel('cubic interpolation error(x)')
	plt.title('error')
	plt.ylim([-1,13])
	plt.grid(True)
	
	plt.subplot(212)
	plt.plot(x, e2, 'ro', linestyle = ":")
	plt.xlabel('x')
	plt.ylabel('lagrange interpolation error(x)')
	plt.ylim([-0.00005,0.00015])
	plt.grid(True)
	plt.show()


def main():
	# print("NCM: Assignment #5: Interpolation \n")
	data = function_values(-3, 4, 15)

	# print("Finding Lagrange interpolation polynom for given data range \n")
	l = lagrange_polynom(data)	
	# print("Lagrange polynom: ", l)

	# print("Finding cubic interpolation polynoms for given data range \n")
	res = get_spline(data)
	# for i in range(len(res)):
	# 	print(res[i], '\n')

	# error arrays
	e1, e2 = interpolation_error(data, res, l)

	# Visualisation	
	plot_interpolate(data, res, l)
	plot_error_interpolation(data, e1, e2)

	
if __name__ == '__main__':
	main()