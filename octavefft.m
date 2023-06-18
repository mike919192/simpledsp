
disp("Octave reading binary data!")

myfile = fopen ("example.bin", "r");
val = fread (myfile, Inf, "double");
fclose(myfile);
data = zeros(1, length(val) / 2);

disp("Octave binary length: ")
disp(length(data))

for i = 1:length(data)
    data(i) = val(((i-1) * 2) + 1) + val(((i-1) * 2) + 2) * i; 
endfor

fft_data = fft(data);

val2 = zeros(1, length(fft_data));

for i = 1:length(fft_data)
    val2(((i-1) * 2) + 1) = real(fft_data(i));
    val2(((i-1) * 2) + 2) = imag(fft_data(i));
endfor

myfile2 = fopen ("example2.bin", "w");
fwrite(myfile2, val2, "double");
fclose(myfile2);

disp("Octave done fft calc!")
