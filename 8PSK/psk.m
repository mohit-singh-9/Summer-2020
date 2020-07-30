

%SNR range
snrlen=10;

%SNR in dB and actual per bit 
%(Check Proakis for factor of 6)
snr_db = linspace(0,snrlen,snrlen);
snr = 6.*10.^(0.1.*snr_db);

%Bitstream size
bitsimlen = 99999;

%Symbol stream size
simlen = bitsimlen / 3;

bits = randi([0,1],1,bitsimlen);


% Converting bits into gray code
sym = bi2de(reshape(bits,3,length(bits)/3).','left-msb');
gray_sym = bitxor(sym,floor(sym/2));
symbol_lst=zeros(simlen,2);
s = zeros(8,2);

for i = 1:8
	s(i,:) = [cos((i-1)*2*pi/8) , sin((i-1)*2*pi/8)];
end
for i = 1:simlen
    symbol_lst(i,:)= s(gray_sym(i)+1,:);
end

symbol_lst=symbol_lst.';
ser=[];
ser_anal=[];
ber=[];

gray = zeros(8,3);
gray(1,:) = [0 0 0];
gray(2,:) = [0 0 1];
gray(3,:) = [0 1 1];
gray(4,:) = [0 1 0];
gray(5,:) = [1 1 0];
gray(6,:) = [1 1 1];
gray(7,:) = [1 0 1];
gray(8,:) = [1 0 0];

  
A= zeros(8,2,2);
A(1,:,:) = [sqrt(2)-1 sqrt(2)-1; 1 -1];
A(2,:,:) = [sqrt(2)+1 -sqrt(2)+1;-1 1];
A(3,:,:) = [-sqrt(2)-1 sqrt(2)+1;1 1];
A(4,:,:) = [sqrt(2)-1 -sqrt(2)-1; 1 -1];
A(5,:,:) = [-sqrt(2)+1 -sqrt(2)+1;-1  1];
A(6,:,:) = [-sqrt(2)-1 sqrt(2)-1;1 -1];
A(7,:,:) = [sqrt(2)+1  -sqrt(2)-1;-1 -1];
A(8,:,:) = [-sqrt(2)+1 sqrt(2)+1;-1  1];

for k = 1:10
    noise = randn(2,simlen);
    y = sqrt(snr(k)).* symbol_lst + noise;
    t=0;
    
    brx=[];
    for i=1:simlen
        z=y(:,i);
        for j=1:8
            y1= [A(j,:,1);A(j,:,2)]*z;
            if and(y1(1)>=0,y1(2)>=0)
                srx = s(j,:);
                break
            end
        end
         brx= [brx,gray(j,:)];
        %for n= 1:8
            %if and(s(n,1)==srx(1), s(n,2) == srx(2))
                %brx= [brx,gray(n,:)];
                %break
            %end
        %end
        
        if and(symbol_lst(1,i)==srx(1),symbol_lst(2,i)==srx(2))
            t=t+1;
        end
    end
  
    ser=[ser,1-(t/simlen)];
    
    ser_anal = [ser_anal,erfc(sqrt(snr(k)/2)*sin(pi/8))];
    %brx = reshape(brx,[1,bitsimlen]);
    bit_diff = bits-brx;
    count = nnz(bit_diff);
    ber = [ber, (count/bitsimlen)]; 
    
    
end

%Plots
close all; figure
semilogy(snr_db,ser_anal,'o--','LineWidth',1)
hold on
semilogy(snr_db,ser,'r--','LineWidth',1)
hold on
semilogy(snr_db,ber,'bs-','LineWidth',1)
hold on
%axis([0 10 0.4 2.5])
grid on
xlabel('SNR in dB')
ylabel('Error Rate')
title('BER & SER curve for 8PSK')





	