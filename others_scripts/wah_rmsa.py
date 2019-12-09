import math.pi
import numpy

#INPUTS: list of pairable vectors. Obj1 should be a list of vectors
#that you want to match to the corresponding obj2.
def vector_angle(v1, v2):
    v1_u = v1 / numpy.linalg.norm(v1)
    v2_u = v2 / numpy.linalg.norm(v2)
    return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))

def zero_vector_pair(v_pair):
    trans = [0. - v_pair[0][0], 0. - v_pair[0][1], 0. - v_pair[0][2]]
    new_pair = []
    for point in v_pair:
        new_point = numpy.array(point) + numpy.array(trans)
        new_pair.append(new_point)
    return new_pair

def calc_rmsa(obj1, obj2, ratio=(pi/6.0)):
    compare_vectors = zip(obj1, obj2)
    vector_ang_sum = 0.0
    for vector_pairs in compare_vectors:
        vector1 = zero_vector_pair(vector_pairs[0])
        vector2 = zero_vector_pair(vector_pairs[1])
        vector_ang_sum += vector_angle(numpy.array(vector1[1]), numpy.array(vector2[1]))
    rmsa = ((vector_ang_sum/ratio)**2 / len(compare_vectors))**0.5
    return rmsa
