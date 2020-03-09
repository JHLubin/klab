# RMSA functions provided by Will Hansen
#INPUTS: list of pairable vectors. Obj1 should be a list of vectors
#that you want to match to the corresponding obj2.

def vector_angle(v1, v2):
	v1_u = v1 / np.linalg.norm(v1)
	v2_u = v2 / np.linalg.norm(v2)
	return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def zero_vector_pair(v_pair):
	trans = [0. - v_pair[0], 0. - v_pair[1], 0. - v_pair[2]]
	new_pair = []
	for point in v_pair:
		new_point = np.array(point) + np.array(trans)
		new_pair.append(new_point)
	return new_pair

def calc_rmsa(obj1, obj2, ratio=(pi/6.0)):
	assert len(obj1) == len(obj2)
	compare_vectors = zip(obj1, obj2)
	vector_ang_sum = 0.0
	for vector_pairs in compare_vectors:
		vector1 = zero_vector_pair(vector_pairs[0])
		vector2 = zero_vector_pair(vector_pairs[1])
		vector_ang_sum += vector_angle(np.array(vector1[1]), np.array(vector2[1]))
	rmsa = ((vector_ang_sum/ratio)**2 / len(obj1))**0.5
	return rmsa


def get_vector_obj_for_rmsa(pose, residue_number):
	"""
	For a given pose and residue number, returns a list of the vectors from CA 
	to N and CA to C.
	"""
	target_res = pose.residue(residue_number)
	CA_N_vector = list(target_res.atom('N').xyz()-target_res.atom('CA').xyz())
	CA_C_vector = list(target_res.atom('C').xyz()-target_res.atom('CA').xyz())

	return [CA_N_vector, CA_C_vector]