"""
This script is a Python conversion of the MATLAB file testing.m.
It includes a proper implementation of the magic function:
 - For n=1: returns a 1x1 magic square.
 - For n=2: raises a ValueError, since a true magic square does not exist.
 - For odd n (n>=3): uses the Siamese method.
 - For doubly even n (n divisible by 4): uses the standard algorithm for doubly even magic squares.

The tester function builds a dictionary with keys 'one', 'two', 'three', and 'four'.
"""

def magic(n):
    """Generate an n x n magic square.

    For n = 1, returns [[1]].
    For n = 2, raises a ValueError since a true magic square does not exist.
    For odd n, uses the Siamese method.
    For doubly even n (n % 4 == 0), uses the standard doubly even algorithm.
    """
    if n == 1:
        return [[1]]
    elif n == 2:
        raise ValueError('Magic square is not possible for n=2')
    elif n % 2 == 1:
        return magic_odd(n)
    elif n % 4 == 0:
        return magic_doubly_even(n)
    else:
        # Singly even numbers (like 6, 10, ...) are not implemented in this version.
        raise NotImplementedError('Magic square generation for singly even numbers is not implemented.')

def magic_odd(n):
    """Generate an odd-order magic square using the Siamese method."""
    # Initialize an n x n square filled with zeros
    square = [[0 for _ in range(n)] for _ in range(n)]
    i, j = 0, n // 2
    for num in range(1, n * n + 1):
        square[i][j] = num
        # Calculate next position
        new_i = (i - 1) % n
        new_j = (j + 1) % n
        if square[new_i][new_j] != 0:
            i = (i + 1) % n
        else:
            i, j = new_i, new_j
    return square

def magic_doubly_even(n):
    """Generate a doubly even (n divisible by 4) magic square."""
    # Fill the square with numbers 1 to n*n sequentially
    square = [[i * n + j + 1 for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            # In a doubly even magic square, if (i % 4 == j % 4) or ((i + j) % 4 == 3), replace the number
            if (i % 4 == j % 4) or ((i + j) % 4 == 3):
                square[i][j] = n * n + 1 - square[i][j]
    return square

def tester(test):
    """Augment the test dictionary with magic squares for orders 1, 2, 3, and 4.
    For n=2, if a magic square is not possible, store None."""
    test['one'] = magic(1)
    try:
        test['two'] = magic(2)
    except ValueError as e:
        test['two'] = None
    test['three'] = magic(3)
    test['four'] = magic(4)
    return test

def testing():
    """Main testing function. Initializes a dictionary, processes it with tester, and prints the results."""
    test = {}
    test = tester(test)

    print('test.one:', test['one'])
    print('test.two:', test['two'])
    print('test.three:')
    if test['three'] is not None:
        for row in test['three']:
            print(' '.join(str(x) for x in row))
    else:
        print(test['three'])

    print('test.four:')
    if test['four'] is not None:
        for row in test['four']:
            print(' '.join(str(x) for x in row))
    else:
        print(test['four'])

if __name__ == '__main__':
    testing()