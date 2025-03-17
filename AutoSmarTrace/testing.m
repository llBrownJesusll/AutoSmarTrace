function testing()

test = struct
test = tester(test)
disp(test.one)
disp(test.two)
disp(test.three)
disp(test.four)
disp(test.five)
end

function test = tester(test)
test.one = magic(1)
test.two = magic(2)
test.three = magic(3)
test.four = magic(4)
test.five = magic(5)
end