function D = eigenproperties_response_luis(ss,n_output,n_input,list_eig)
%% inputs:
%ss: state-space model
%n_output: observabibility in this output
%n_input: observabibility in this input
%list_eig: list of number of eigenvalues considered e.g. [1:1:44] or [1:1:10
%12:1:44]--> so in this last case the eigenvalue 11 is not considered
%% outputs:
% figure with the response
%% Code:
%D=eig(ss.A)
[R,D,L] = eig(ss.A);
D=diag(D);
B = inv(R)*ss.B;
C = ss.C*R;
D_ss = ss.D(n_output,n_input);

i=1;

t_final = 1;
t_step  = 1e-4;

time = 0:t_step:t_final;
response = zeros(size(time));

for t = time %time response
    response(i)=0;
    for mode=list_eig
        eigenvalue = D(mode);         
        residue = C(n_output,mode)*B(mode,n_input);
        response(i) = response(i) + real( (residue/eigenvalue) * ( exp(eigenvalue*t) -1 ) ) ;
    end
    response(i) = response(i) + D_ss;
    i=i+1;
end
figure
plot([0:t_step:t_final],response)
hold on
y=step(ss,[0:t_step:t_final]);
plot([0:t_step:t_final],y(:,n_output,n_input),'Color','red','LineStyle','--')
grid on
legend('Decomposed response','State-space response')
xlabel('Time [s]')
ylabel('Output')
end



