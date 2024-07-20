"""
tests/test_jacobian
"""

import pprint
from besmarts.core import geometry

# x1 = [[0,1,0]], [[0,0,0]], [[1,0,0]]
x2 = [[.5,1,0]], [[.5,0,0]], [[1,.2,-.1]]

result = geometry.jacobian2_angle(*x2)[0]

for row in geometry.jacobian2_reshape_to_matrix(result):
    print(*map("{: 11.4e}".format, row))

"""
 1.5086e-02  9.8058e-01  7.5429e-02  6.0343e-02 -9.8058e-01  3.0172e-01 -7.5429e-02  2.2204e-16 -3.7715e-01
 9.8058e-01 -5.9629e-17 -1.9612e-01 -9.8058e-01  5.9629e-17  1.9612e-01  0.0000e+00  0.0000e+00  0.0000e+00
 7.5429e-02 -1.9612e-01  3.7715e-01  3.0172e-01  1.9612e-01  1.5086e+00 -3.7715e-01 -2.7756e-17 -1.8857e+00
 6.0343e-02 -9.8058e-01  3.0172e-01 -2.2646e+00  3.3776e+00  8.3810e-03  2.2042e+00 -2.3970e+00 -3.1010e-01
-9.8058e-01  5.9629e-17  1.9612e-01  3.3776e+00  2.2662e+00 -6.7551e-01 -2.3970e+00 -2.2662e+00  4.7939e-01
 3.0172e-01  1.9612e-01  1.5086e+00  8.3810e-03 -6.7551e-01 -2.2243e+00 -3.1010e-01  4.7939e-01  7.1574e-01
-7.5429e-02  0.0000e+00 -3.7715e-01  2.2042e+00 -2.3970e+00 -3.1010e-01 -2.1288e+00  2.3970e+00  6.8724e-01
 2.2204e-16  0.0000e+00 -2.7756e-17 -2.3970e+00 -2.2662e+00  4.7939e-01  2.3970e+00  2.2662e+00 -4.7939e-01
-3.7715e-01  0.0000e+00 -1.8857e+00 -3.1010e-01  4.7939e-01  7.1574e-01  6.8724e-01 -4.7939e-01  1.1700e+00
"""