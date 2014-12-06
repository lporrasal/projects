def find_vowels("aaayhuiommmm"):

    count = 0
    vowels = "aeiuoAEIOU"
    for letter in sentence:
        if letter in vowels:
            count += 1
    print count

if __name__ == '__main__':
    import doctest
    doctest.testmod("aaayhuiommmm")