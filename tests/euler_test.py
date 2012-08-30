from quaternion import Quaternion
import itertools
import math

orders = [
    'xyz', 'xzy', 'xyx', 'xzx',
    'yzx', 'yxz', 'yxy', 'yzy',
    'zxy', 'zyx', 'zxz', 'zyz']

test_angles = range(-360, 360, 45)

def check_euler(order, angles):
    q = Quaternion.from_euler(order, map(math.radians, angles))
    euler = q.to_euler(order)
    r = Quaternion.from_euler(order, euler)
    if (r.dot(q) < 0):
        r.negate()
    for a, b in zip(q.components, r.components):
        assert abs(a - b) < 1e-6, "%s != %s" % (str(r), str(q))
    
def test_euler():
    for angles in itertools.combinations(test_angles, 3):
        for order in orders:
            yield check_euler, order, angles
