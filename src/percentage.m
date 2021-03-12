indices = [1:15962 (29525+5000):72540 (86242+5000):129072 (142960+5000):186050 (199677+5000):242583 (256394+5000):299561 (313112+5000):332701];
error = instant_error_perform(indices,:);


x = abs(error(:,1)) <= 20;
y = abs(error(:,2)) <= 20;
z = abs(error(:,3)) <= 20;
total1 = 0;
total2 = 0;
total3 = 0;
total = 0;
total4 = 0;
for i=1:length(y)
    total1 = total1 + x(i);
    total2 = total2 + y(i);
    total3 = total3 + z(i);
    total = total + y(i)*z(i);
    total4 = total4 + x(i)*y(i)*z(i);
end

per1 = total1/length(y) %X < 20
per2 = total2/length(y) %Y < 20
per3 = total3/length(y) % Z < 20
per = total/length(y) %YZ < 20
per4 = total4/length(y) %XYZ < 20
