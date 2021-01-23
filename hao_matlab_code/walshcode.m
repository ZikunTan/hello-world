function [op_code] = walshcode(a)
len_walsh = 2 .^ ceil(log2(a));  %length of walsh code , and number of walsh codes
walsh_code = 1;
i = 2;
while (i <= len_walsh)
    walsh_code = [walsh_code, walsh_code ; walsh_code, -1.* walsh_code];
    i = i.*2;
end
op_code = walsh_code;
%if you want to get only one sequence then remove the '%' from the below
%line and put it in the begining of the above
%this fucntion randomly picks one row(one sequence) from the code matrix
%op_code = ranpick(len_walsh, walsh_code);
% function [op_walsh] = ranpick(a,b)
% code_walsh = b;
% op_walsh = [];
% ran_pick = 10 .* rand;
% ran_pick = floor(ran_pick + 1)
% no_seq = a;
% if ran_pick <= no_seq
%     code_walsh
%     op_walsh = code_walsh(ran_pick, :)
% else
%     op_walsh = ranpick(no_seq, code_walsh);
% end