indices = [1:22092 (43496+15000):78813 (100217+15000):135533 (156937+15000):192189 (213657+15000):248589 (270377+15000):305525 (327097+15000):332701];
error = instant_error_perform(:,:);


x = abs(error(:,1)) <= 20;
y = abs(error(:,2)) <= 20;
total1 = 0;
total2 = 0;
total = 0;
for i=1:length(y)
    total1 = total1 + x(i);
    total2 = total2 + y(i);
    total = total + y(i)*x(i);
end

per1 = total1/length(y) %X < 20
per2 = total2/length(y) %Y < 20
per = total/length(y) %YZ < 20
