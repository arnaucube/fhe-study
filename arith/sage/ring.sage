def run_test(a, b):
	print("\nnew test:")
	print(a)
	print(b)
	c = a*b
	print(c)
	print(c.list())



n_iters = 100
Q= 65537
print(Q)
N=4
F = GF(Q)
R = QuotientRing(F[x], x^N + 1, names="X")
print(R)

a = R([4,2,1,0])
b = R([1,2,3,4])
run_test(a,b)


#  print("Elements of the polynomial ring:")
#  for e in R:
#      print(e)


# Other:
# ======
#
#  t = R.gen()
#  a =  0 + t   + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5
#  b =  5 + 4*t + 3*t^2 + 2*t^3 + 1*t^4 + 0*t^5
#  print("add", a+b)
#  print("sub", a-b)
#  print("mul", a*b)
#  a =  0 + t   + 2*t^2 + 3*t^3 + 4*t^4 + 5*t^5

# print("ring elem mul testvectors")
#
#
# def randvec(size=N):return [int(random()*(Q-1)) for t in range(size)]
#
# a_vecs = [None]*n_iters
# b_vecs = [None]*n_iters
# c_vecs = [None]*n_iters
#
# for i in range(n_iters):
#     a_vec = randvec()
#     b_vec = randvec()
#     a_pol = R(a_vec)
#     b_pol = R(b_vec)
#
#     c_pol = a_pol*b_pol
#
#     a_vecs[i] = a_pol.list()
#     b_vecs[i] = b_pol.list()
#     c_vecs[i] = c_pol.list()
#
# print("let a_vecs = vec!{};\n".format(a_vecs))
# print("let b_vecs = vec!{};\n".format(b_vecs))
# print("let c_vecs = vec!{};".format(c_vecs))


# # cyclotomic
# 
# Q= 65537
# print(Q)
# N=4
# F = GF(Q)
# R = QuotientRing(F[x], x^N - 1, names="X")
# print(R)
# 
# a = R([1,0,0,2])
# b = R([0,0,0,2])
# run_test(a, b)
