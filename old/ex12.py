print "How old are you?",
age = int(raw_input())
# one can also use age = raw_input("How old are you? ")
print "How tall are you?",
height = int(raw_input())
print "How much do you weigh?",
weight = raw_input()

print "So, you're %r old, %r tall and %r heavy." % (
    age, height, weight)

print age + height