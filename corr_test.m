% x=1:100;
% % y=(5+1.3*x).*(rand(size(x))*0.2+0.5);
% y=(10+1.3*x).*(1+randn(size(x))*0.2);

n=200;
x=[1:n];
x=x.*(rand(size(x))*0.2 +0.5);
y=[1:n];
y=y.*(rand(size(y))*0.2 +0.5);



% Fit line to data.
[fitresult, gof] = fit( x', y', 'poly1' )

% Plot fit with data.
figure( 'Name', 'Line fit' );
h = plot( fitresult, x', y' );
legend( h, 'data', 'linear fit', 'Location', 'NorthEast' );
% Label axes
xlabel x
ylabel y
grid on

