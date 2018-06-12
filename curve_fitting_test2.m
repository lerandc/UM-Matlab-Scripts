clf
clearvars
func = @(x,A,B,C,D,E) D.*(x-E).^2+A.*(x-B)+C;
iter_sum = 0;
r_sum = 0;

    clf 
    a = randi([-100 100],1);
    b = randi([-100 100],1);
    c = randi([-100 100],1);
    d = 0;%randi([-100 100],1);
    e = randi([-100 100],1);
    
    x_test = linspace(-100,100,1000);
    y_test = func(x_test,a,b,c,d,e);
    plot(x_test,y_test);
    grid on

    ga = a + 1*randn(1);
    gb = b + 1*randn(1);
    gc = c + 1*randn(1);
    gd = a + 1*randn(1);
    ge = b + 1*randn(1);
    
    y_test = y_test + (5.*10.^((log10(abs(d))))).*((randn([1,1000])));
    tic
    [A,B,C,D,E,iter,residual] = fit_algorithm(ga,gb,gc,gd,ge,x_test,y_test)
    toc
    
function [fitA, fitB, fitC,fitD,fitE,iter,residual] = fit_algorithm(ga,gb,gc,gd,ge,x0,y0)
    func = @(x,A,B,C,D,E) D.*(x-E).^2+A.*(x-B)+C;
    y_fit = func(x0,ga,gb,gc,gd,ge);
    residual = sum((abs(y_fit - y0).^2));
    tolerance = 1e-2;
    dV = 10;
    A = ga;
    B = gb; 
    C = gc;
    D = gd;
    E = ge;
    iter = 0;
    r_track = 10;
    r_hold = residual;
    while((r_track > tolerance) && (iter < 100000) && (dV > 1e-4))
        r1 = sum((abs(func(x0,A-dV,B,C,D,E)-y0).^2));
        r2 = sum((abs(func(x0,A+dV,B,C,D,E)-y0).^2));
        r3 = sum((abs(func(x0,A,B-dV,C,D,E)-y0).^2));
        r4 = sum((abs(func(x0,A,B+dV,C,D,E)-y0).^2));
        r5 = sum((abs(func(x0,A,B,C-dV,D,E)-y0).^2));
        r6 = sum((abs(func(x0,A,B,C+dV,D,E)-y0).^2));
        r7 = sum((abs(func(x0,A,B,C,D-dV,E)-y0).^2));
        r8 = sum((abs(func(x0,A,B,C,D+dV,E)-y0).^2));
        r9 = sum((abs(func(x0,A,B,C,D,E-dV)-y0).^2));
        r10 = sum((abs(func(x0,A,B,C,D,E+dV)-y0).^2));
        [R,I] = min([r1;r2;r3;r4;r5;r6;r7;r8;r9;r10]);
        
        if R > residual
            dV = dV/1.5;
        else
            residual = R;
            switch I
                case 1
                    A = A-dV;
                case 2
                    A = A+dV;
                case 3
                    B = B-dV;
                case 4
                    B = B+dV;
                case 5
                    C = C-dV;
                case 6
                    C = C+dV;
                case 7
                    D = D-dV;
                case 8
                    D = D+dV;
                case 9
                    E = E-dV;
                case 10
                    E = E+dV;
            end
           clf
            plot(x0,y0,'.','MarkerSize',3)
            hold on
            plot(x0,func(x0,A,B,C,D,E),'--')
            grid on
            drawnow
            iter = iter + 1;
            r_track = r_hold - residual;
            r_hold = residual;
        end
    end
    fitA = A;
    fitB = B; 
    fitC = C;
    fitD = D;
    fitE = E;
    plot(x0,y0,'.','MarkerSize',3)
    hold on
    plot(x0,func(x0,A,B,C,D,E),'--R')
    grid on
end
