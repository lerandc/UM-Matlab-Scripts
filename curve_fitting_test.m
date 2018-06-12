clf
clearvars
func = @(x,A,B,C) A.*(x-B)+C;
iter_sum = 0;
r_sum = 0;

    clf 
    a = randi([-100 100],1);
    b = randi([-100 100],1);
    c = randi([-100 100],1);

    x_test = linspace(-100,100,1000);
    y_test = func(x_test,a,b,c);
    plot(x_test,y_test);
    grid on

    ga = a + 100*randn(1);
    gb = b + 100*randn(1);
    gc = c + 100*randn(1);
    y_test = y_test + (1.*10.^((log10(abs(a))))).*((randn([1,1000])));
    tic
    [A,B,C,iter,residual] = fit_algorithm(ga,gb,gc,x_test,y_test)
    toc
    plot(x_test,func(x_test,a,b,c),'--K')
function [fitA, fitB, fitC,iter,residual] = fit_algorithm(ga,gb,gc,x0,y0)
    func = @(x,A,B,C) A.*(x-B)+C;
    y_fit = func(x0,ga,gb,gc);
    residual = sum((abs(y_fit - y0).^2));
    tolerance = 1e-3;
    dV = 10;
    A = ga;
    B = gb; 
    C = gc;
    iter = 0;
    r_track = 1;
    r_hold = residual;
    while((r_track > tolerance) && (iter < 10000) && (dV > 1e-4))
        r1 = sum((abs(func(x0,A-dV,B,C)-y0).^2));
        r2 = sum((abs(func(x0,A+dV,B,C)-y0).^2));
        r3 = sum((abs(func(x0,A,B-dV,C)-y0).^2));
        r4 = sum((abs(func(x0,A,B+dV,C)-y0).^2));
        r5 = sum((abs(func(x0,A,B,C-dV)-y0).^2));
        r6 = sum((abs(func(x0,A,B,C+dV)-y0).^2));
        [R,I] = min([r1;r2;r3;r4;r5;r6]);
        if R > residual
            dV = dV/2;
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
            end
           clf
            plot(x0,y0,'.','MarkerSize',3)
            hold on
            plot(x0,func(x0,A,B,C),'--')
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
    plot(x0,y0,'.','MarkerSize',3)
    hold on
    plot(x0,func(x0,A,B,C),'R')
    grid on
end
