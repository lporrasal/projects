class Animal:
	def __init__(self, name, age):
		self.name = name
		self.age = age

	def bark(self):
		print("My dog {0} is {1} years old, they say {2}!\n".format(self.name, self.age, "guau"))

dog = Animal("Rarra",7)
dog.bark()