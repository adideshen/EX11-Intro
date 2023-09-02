from itertools import combinations
import copy

class Node:
	def __init__(self, data, positive_child=None, negative_child=None):
		self.data = data
		self.positive_child = positive_child
		self.negative_child = negative_child


class Record:
	def __init__(self, illness, symptoms):
		self.illness = illness
		self.symptoms = symptoms
	
			
def parse_data(filepath):
	with open(filepath) as data_file:
		records = []
		for line in data_file:
			words = line.strip().split()
			records.append(Record(words[0], words[1:]))
		return records
		
		
class Diagnoser:
	def __init__(self, root: Node):
		self.root = root

	def diagnose(self, symptoms: list):
		"""
		The function diagnoses an illness based on a list of symptoms
		:param symptoms: List of symptoms
		:return: The diagnosed illness
		"""

		if self.root.negative_child is None:
			return self.root.data
		else:
			if self.root.data in symptoms:
				positive_header = Diagnoser(self.root.positive_child)
				return positive_header.diagnose(symptoms)
			else:
				negative_header = Diagnoser(self.root.negative_child)
				return negative_header.diagnose(symptoms)

	def calculate_success_rate(self, records) -> float:
		"""
		The function calculates how many diseases it has been able to diagnose from a list of diseases and symptoms.
		:param records: List of objects from the class Record.
		:return: The ratio between the number of successes of the tree on the records in the records and the total
		number of records
		"""
		num_of_successes = 0
		for i in range(len(records)):
			if records[i].illness == self.diagnose(records[i].symptoms):
				num_of_successes += 1
		if len(records) == 0:
			raise ValueError("There are no symptoms")
		else:
			return num_of_successes/len(records)

	def helper_all_illnesses(self, dict_of_trees_illnesses):
		"""
		The function calculates the number of times each of the trees illnesses appear in the tree
		:param dict_of_trees_illnesses: Empty dictionary
		:return: A dictionary that its values are the illnesses that appear in the tree and the number of times
		they appear
		"""
		if not self.root.negative_child:
			if self.root.data in dict_of_trees_illnesses:
				dict_of_trees_illnesses[self.root.data] += 1
			else:
				dict_of_trees_illnesses[self.root.data] = 1
		else:
			negative_diagnoser = Diagnoser(self.root.negative_child)
			positive_diagnoser = Diagnoser(self.root.positive_child)
			negative_diagnoser.helper_all_illnesses(dict_of_trees_illnesses)
			positive_diagnoser.helper_all_illnesses(dict_of_trees_illnesses)
		return dict_of_trees_illnesses

	def all_illnesses(self):
		"""
		:return: A list of all the illnesses that appear in the tree which are arranged according to the number of
		times they appear
		"""
		dict_illnesses = {}
		full_dict_illnesses = self.helper_all_illnesses(dict_illnesses)
		list_of_tuples_of_illnesses = sorted(full_dict_illnesses.items(), key=lambda x: x[1], reverse=True)
		list_of_illnesses = []
		for illnesses in list_of_tuples_of_illnesses:
			list_of_illnesses.append(illnesses[0])
		return list_of_illnesses

	def helper_paths_to_illness(self, illness, current_path, all_paths):
		"""
		This function find all the paths to an illness in a tree
		:param illness: The illness that its paths need to be found.
		:param current_path: A single path that is found.
		:param all_paths: List of all of the paths.
		:return: The list of all of the paths.
		"""
		if self.root.negative_child is None:
			if self.root.data == illness:
				all_paths.append(current_path)
		else:
			positive_child = Diagnoser(self.root.positive_child)
			negative_child = Diagnoser(self.root.negative_child)
			return negative_child.helper_paths_to_illness(illness, current_path + [False], all_paths), positive_child.helper_paths_to_illness(illness, current_path + [True], all_paths)

	def paths_to_illness(self, illness):
		"""
		This function find all the paths to an illness in a tree
		:param illness: The illness that its paths need to be found.
		:return: List of paths in a tree.
		"""
		all_paths = []
		self.helper_paths_to_illness(illness,[], all_paths)
		return all_paths


	def are_my_children_leaves(self):
		"""
		This function checks if ones children are leaves.
		:return: True is they are False if they not
		"""
		if self.root.negative_child:
			positive_child = self.root.positive_child
			negative_child = self.root.negative_child
			if positive_child.positive_child is None and negative_child.positive_child is None:
				return True
			else:
				return False
		else:
			return False

	def minimize_for_false(self):
		"""
		This function minimize a tree if one of the vertices has two sons with the same data
		:return: None
		"""
		if self.are_my_children_leaves():
			if self.root.negative_child is not None and self.root.positive_child is not None:
				if self.root.negative_child.data == self.root.positive_child.data:
					self.root.data = self.root.positive_child.data
					self.root.positive_child = None
					self.root.negative_child = None
				else:
					return
		elif self.root.negative_child is None:
			return
		else:
			positive_child = Diagnoser(self.root.positive_child)
			negative_child = Diagnoser(self.root.negative_child)
			positive_child.minimize_for_false()
			negative_child.minimize_for_false()

	def minimize_for_true(self):
		"""
		This function minimize a tree if one of the vertices has no children
		:return: None
		"""
		if self.root.negative_child.data is None:
			self.root.data = self.root.positive_child.data
			self.root.negative_child = self.root.positive_child.negative_child
			self.root.positive_child = self.root.positive_child.positive_child
		elif self.root.positive_child.data is None:
			self.root.data = self.root.negative_child.data
			self.root.positive_child = self.root.negative_child.positive_child
			self.root.negative_child = self.root.negative_child.negative_child
		elif self.are_my_children_leaves():
			if self.root.positive_child.data is None:
				self.root.data = self.root.negative_child.data
				self.root.positive_child = self.root.negative_child.positive_child
				self.root.negative_child = self.root.negative_child.negative_child
			elif self.root.negative_child.data is None:
				self.root.data = self.root.positive_child.data
				self.root.negative_child = self.root.positive_child.negative_child
				self.root.positive_child = self.root.positive_child.positive_child
		else:
			positive_child = Diagnoser(self.root.positive_child)
			negative_child = Diagnoser(self.root.negative_child)
			positive_child.minimize_for_true()
			negative_child.minimize_for_true()

	def minimize(self, remove_empty=False):
		"""
		This function minimize a tree if one of its vertices is unnecessary.
		:param remove_empty: True if it needs to reduce the leaves of a vertices that both of it vertices are None.
		False if not.
		:return: None
		"""
		self.minimize_for_false()
		if remove_empty:
			self.minimize_for_true()
		return


def appropriate_record(records, positive_symptoms, negative_symptoms):
	"""
	The function receives conditions for finding a record to place in a leaf returns the most appropriate illness
	:param records: List of all the records.
	:param positive_symptoms: The symptoms appear in the leaf's illness
	:param negative_symptoms: The symptoms that do not appear in the leaf's illness
	:return: The illness that is suitable to be placed in the leaf.
	"""
	all_appropriate_illnesses = {}
	for i in range(len(records)):
		record_symptoms = records[i].symptoms
		fine_record = True
		for j in range(len(positive_symptoms)):
			if positive_symptoms[j] not in record_symptoms:
				fine_record = False
				break
		if fine_record:
			for k in range(len(negative_symptoms)):
				if negative_symptoms[k] in record_symptoms:
					fine_record = False
					break
		if fine_record:
			if records[i].illness in all_appropriate_illnesses:
				all_appropriate_illnesses[records[i].illness] += 1
			else:
				all_appropriate_illnesses[records[i].illness] = 1
	dict(sorted(all_appropriate_illnesses.items(), key=lambda item: item[1]))
	if len(all_appropriate_illnesses) == 0:
		return None
	most_common_illness = max(all_appropriate_illnesses, key=all_appropriate_illnesses.get)
	return most_common_illness


def my_tree(records, symptoms, node, counter, positive_symptoms, negative_symptoms):
	"""
	This function builds a tree.
	:param records: List of records
	:param symptoms: A list of symptoms that will appear in the tree
	:param node: The tree root
	:param counter:  The location of the symptom in the list to be embedded in the tree.
	:param positive_symptoms: List of the symptoms appear in the leaf's illness
	:param negative_symptoms: List of the symptoms that do not appear in the leaf's illness
	:return: None
	"""
	if counter == len(symptoms):
		illness = appropriate_record(records, positive_symptoms, negative_symptoms)
		if positive_symptoms and node.data == positive_symptoms[-1]:
			node.positive_child = Node(illness)
			symptom = positive_symptoms.pop()
			negative_symptoms.append(symptom)
			illness2 = appropriate_record(records, positive_symptoms, negative_symptoms)
			node.negative_child = Node(illness2)
		elif negative_symptoms and node.data == negative_symptoms[-1]:
			node.negative_child = Node(illness)
			symptom =  negative_symptoms.pop()
			positive_symptoms.append(symptom)
			illness2 = appropriate_record(records, positive_symptoms, negative_symptoms)
			node.positive_child = Node(illness2)
	else:
		current_symptom = symptoms[counter]
		node.positive_child = Node(current_symptom)
		node.negative_child = Node(current_symptom)
		positive_symptoms.append(symptoms[0])
		my_tree(records, symptoms, node.positive_child, counter + 1, (positive_symptoms + [current_symptom])[:] , negative_symptoms[:])
		negative_symptoms.append(positive_symptoms.pop())
		my_tree(records, symptoms, node.negative_child, counter + 1, positive_symptoms[:], (negative_symptoms + [current_symptom])[:])


def most_comment_illness(records):
	"""
	:param records: List of records
	:return: The illness that is most often on the list
	"""
	all_illnesses = {}
	for i in range(len(records)):
		if records[i].illness in all_illnesses:
			all_illnesses[records[i].illness] += 1
		else:
			all_illnesses[records[i].illness] = 1
	dict(sorted(all_illnesses.items(), key=lambda item: item[1]))
	most_common_illness = max(all_illnesses, key=all_illnesses.get)
	return most_common_illness


def build_tree(records, symptoms):
	"""
	This function builds a tree and returns exceptions if they exist
	:param records: List of records.
	:param symptoms: List of symptoms
	:return: The root of the tree
	"""
	if len(symptoms) == 0:
		diagnoser = Diagnoser(Node(most_comment_illness(records)))
		return diagnoser
	if len(records) == 0:
		diagnoser = Diagnoser(Node(None))
		return diagnoser
	for i in range(len(records)):
		if type(records[i]) != Record:
			raise TypeError("One of the records is not from type Record")
	for i in range(len(symptoms)):
		if type(symptoms[i]) != str:
			raise TypeError("One of the symptoms is not from type str")
	if len(symptoms) > 0:
		if len(symptoms) == 1:
			diagnoser = Diagnoser(Node(symptoms[0]))
			diagnoser.root.negative_child = Node(appropriate_record(records, [], [symptoms[0]]))
			diagnoser.root.positive_child = Node(appropriate_record(records, [symptoms[0]], []))
		else:
			diagnoser = Diagnoser(Node(symptoms[0]))
			my_tree(records, symptoms, diagnoser.root, 1, [], [])
		return diagnoser

def optimal_tree(records, symptoms, depth):
	"""
	This function checks which tree from the a given depth sized tree collection is the tree with the highest success rate
	:param records: List of objects from type Record
	:param symptoms: List of symptoms
	:param depth: The desired tree depth
	:return: The root of a tree from type Diagnose
	"""
	if depth == 0:
		return build_tree(records, [])
	elif depth < 0 or depth > len(symptoms):
		raise ValueError("Depth is not valid")
	for symptom in symptoms:
		if symptoms.count(symptom) > 1:
			raise ValueError("The symptoms list contains double values")
	for i in range(len(records)):
		if type(records[i]) != Record:
			raise TypeError("One of the records is not from type Record")
	for i in range(len(symptoms)):
		if type(symptoms[i]) != str:
			raise TypeError("One of the symptoms is not from type str")
	best_success = (0, 0)
	symptoms_lists = []
	combinations_list = list(combinations(symptoms, depth))
	for x in combinations_list:
		symptoms_lists.append(x)
		tree_root = build_tree(records, x)
		tree_success = tree_root.calculate_success_rate(records)
		if tree_success > best_success[1]:
			best_success = (tree_root, tree_success)
	return best_success[0]


if __name__ == "__main__":
	# covid_leaf1 = Node('covid-19')
	# covid_leaf2 = Node('covid-19')
	# flu_leaf = Node("insomnia", covid_leaf1, covid_leaf2)
	# cold_leaf = Node("cold")
	# fever = Node("fever", flu_leaf, cold_leaf)
	# none_node = Node(None)
	# healthy_leaf = Node('healthy')
	# sweat_leaf = Node("sweat", none_node, healthy_leaf)
	# tree_root = Node("cough", fever, sweat_leaf)
	a = Node("a")
	z = Node("z")
	y = Node("y", z, a)
	none = Node(None)
	x = Node("x", none, y)
	diagnose = Diagnoser(x)
	diagnose.minimize(remove_empty=True)
	x = 0
	# record1 = Record('influenza', ['cough', 'fatigue', 'headache', 'nausea'])
	# record2 = Record('influenza',['headache', 'fever', 'irritability', 'rigidity'])
	# record3 = Record("healthy", [])
	# record4 = Record('cold', ['cough', 'fever'])
	# records = [record3, record2, record1, record4]
	# symptoms = ["headache", "fever"]
	# diagnoser = build_tree(records, symptoms)
	# diagnoser.minimize(remove_empty=True)

	# Manually build a simple tree.
	#                cough
	#          Yes /       \ No
	#        fever           healthy
	#   Yes /     \ No
	# covid-19   cold

	# flu_leaf = Node("covid-19", None, None)
	# cold_leaf = Node("cold", None, None)
	# inner_vertex = Node("fever", flu_leaf, cold_leaf)
	# healthy_leaf = Node("healthy", None, None)
	# root = Node("cough", inner_vertex, healthy_leaf)
	#
	# diagnoser = Diagnoser(root)

	# # Simple test
	# diagnosis = diagnoser.diagnose(["cold"])
	# if diagnosis == "cold":
	# 	print("Test passed")
	# else:
	# 	print("Test failed. Should have printed cold, printed: ", diagnosis)
	#
	# pathes = diagnoser.all_illnesses()
	# print(pathes)
	#
	# Add more tests for sections 2-7 here.