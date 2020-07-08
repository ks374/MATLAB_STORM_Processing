%%
base_path = 'Z:\Chenghang\7.7.20.D2OTesting\D2O\';
outpath = [base_path 'Fitting_result\'];

name_channel = '488';
name_number = '0002';
name = [name_channel 'storm_' name_number '.txt'];
outname = [name_channel 'storm_' name_number '.mat'];
A = textread([outpath name]);
frame_num = A(1,1);
A(1,:) = [];
A = sortrows(A,4);
init_frame = A(1,4);
if name_channel == '647'
    freq = 100;
else
    freq = 60;
end
time_length = floor(frame_num/60) + 1;
%
Height = zeros(time_length,1);
Sum = zeros(time_length,1);
Molecule = zeros(time_length,1);
%
for i = init_frame:freq:(frame_num+frame_num)
    Background_temp = 0;
    Height_temp = 0;
    Sum_temp = 0;
    Molecule_temp = 0;
    while A(1,4) <= (i+59)
        Background_temp = Background_temp + A(1,1);
        Height_temp = Height_temp + A(1,2);
        Sum_temp = Sum_temp + A(1,3);
        Molecule_temp = Molecule_temp + 1;
        A(1,:) = [];
        if isempty(A) == 1
            break
        end
    end
    a = (i - init_frame) / freq + 1;
    disp(a)
    Height(a) = Height_temp / Background_temp;
    Sum(a) = Sum_temp / Background_temp;
    Molecule(a) = Molecule_temp;
    if isempty(A) == 1
        break
    end
end
save([outpath outname],'Height','Sum','Molecule');