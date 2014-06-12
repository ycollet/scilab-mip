A = [1 2 3; 4 5 6; 7 8 9];

serialize_set(A);

B = serialize_get();

disp(A)
disp(B)
