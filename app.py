from flask import Flask, render_template, request
import numpy as np
import sympy
from sympy import *
import math

app = Flask(__name__)
# Define the function whose roots are required
def f(val,func):
    x=val
    return eval(func)


####  NEWTON RAPHSON METHOD  ####
newton_raphs_data = {'iter':[] , 'x':[] , 'f(x)':[] , "f'(x)":[]}
def newton_raphs_f(val,func):
    x = val
    return eval(func)

def derivative(expr):
    x, y = symbols('x y')
    # Use sympy.Derivative() method
    expr_diff = Derivative(expr, x)
    return format(expr_diff.doit())

def newtonRaphson(func,x0,e,N):
    step = 1
    flag = 1
    condition = True
    while condition:
        if newton_raphs_f(x0,derivative(func)) == 0.0:
            error = 'Error! Division By Zero.'
            print(error)
            print(error)
            return render_template('index.html' , error=error)
        
        x1 = x0 - newton_raphs_f(x0,func)/newton_raphs_f(x0,derivative(func))
        print('Iteration-%d, x1 = %0.6f , f(x1) = %0.6f and f\'(x) = %0.6f' % (step, x1, newton_raphs_f(x1,func) , newton_raphs_f(x0,derivative(func))))
        newton_raphs_data['iter'].append(step)
        newton_raphs_data['x'].append(round(x1,4))
        newton_raphs_data['f(x)'].append(round(newton_raphs_f(x1,func),4))
        newton_raphs_data["f'(x)"].append(round(newton_raphs_f(x0,derivative(func)),4))
        x0 = x1
        step = step + 1
        
        if step > N:
            flag = 0
            break
        
        condition = abs(newton_raphs_f(x1,func)) > e
    
    if flag==1:
        print('\nRequired root is: %0.8f' % x1)
        return x1
    else:
        print('\nNot Convergent.')
################  Newton Raphson Ends  #################

################ Secant Methods Functions  #############
secant_data = {'iter':[],'x2':[],'f(x2)':[]}
def secant(x0,x1,e,N,expr):
    print('\n\n*** SECANT METHOD IMPLEMENTATION ***')
    step = 1
    condition = True
    while condition:
        if f(x0,expr) == f(x1,expr):
            print('Divide by zero error!')
            break
        
        x2 = x0 - (x1-x0)*f(x0,expr)/( f(x1,expr) - f(x0,expr) )

        #print('Iteration-%d, x2 = %0.6f and f(x2) = %0.6f' % (step, x2, f(x2,expr)))
        secant_data['iter'].append(step)
        secant_data['x2'].append(round(x2,4))
        secant_data['f(x2)'].append(round(f(x2,expr),4))
        x0 = x1
        x1 = x2
        step = step + 1
        
        if step > N:
            print('Not Convergent!')
            break
        
        condition = abs(f(x2,expr)) > e
    return x2
###########   Secant Method Ends    #############

##########   Regula Falsi Method    #############
falsi_data = {'iter':[],'x2':[],'f(x2)':[]}
def falsi(x0,x1,e,expr):
    step = 1
    condition = True
    while condition:
        x2 = x0 - (x1-x0) * f(x0,expr)/( f(x1,expr) - f(x0,expr) )
        #print('Iteration-%d, x2 = %0.6f and f(x2) = %0.6f' % (step, x2, f(x2,expr)))
        falsi_data['iter'].append(step)
        falsi_data['x2'].append(round(x2,4))
        falsi_data['f(x2)'].append(round(f(x2,expr),4))

        if f(x0,expr) * f(x2,expr) < 0:
            x1 = x2
        else:
            x0 = x2

        step = step + 1
        condition = abs(f(x2,expr)) > e

    print('\nRequired root is: %0.8f' % x2)
    return x2
##########   Regula Falsi Method Ends   ##############
# Re-writing f(x)=0 to x = g(x)
def g(x):
    return 1/math.sqrt(1+x)

############  Fixedpoint Simple Iterartive  ############
fixedpoint_data = {'iter':[] , 'x':[] , 'f(x)':[]}
# Implementing Fixed Point Iteration Method
def fixedpoint(expr,x0, e, N):
    step = 1
    flag = 1
    condition = True
    while condition:
        x1 = g(x0)
        print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1,expr)))
        fixedpoint_data['iter'].append(step)
        fixedpoint_data['x'].append(round(x1,4))
        fixedpoint_data['f(x)'].append(round(f(x1,expr),4))
        x0 = x1

        step = step + 1
        
        if step > N:
            flag=0
            break
        
        condition = abs(f(x1,expr)) > e

    if flag==1:
        return x1
    else:
        return "not convergent"

##########  Fixedpoint Simple Iterartive Method Ends ############

##########  Euler Method Starts Here  #############
euler_data = {'iter':[] , 'x0':[] , 'y0':[] , 'slope':[] , 'yn':[]}
# function to be solved
def euler_f(func,val_x,val_y):
    x = val_x
    y = val_y
    return eval(func)

## Actual EUler Function and Logic
def euler(func,x0,y0,xn,n,input_data):
    
    # Calculating step size
    h = (xn-x0)/n
    
    print('\n-----------SOLUTION-----------')
    print('------------------------------')    
    print('x0\ty0\tslope\tyn')
    print('------------------------------')

    # initial_x0 = x0
    # initial_y0 = y0
    # function = func
    # iter = n
    # calc_point = xn

    for i in range(n):
        slope = euler_f(func,x0, y0)
        yn = y0 + h * slope
        print('%d\t%.4f\t%.4f\t%0.4f\t%.4f'% (i+1,x0,y0,slope,yn) )
        euler_data['iter'].append(i+1)
        euler_data['x0'].append(round(x0,4))
        euler_data['y0'].append(round(y0,4))
        euler_data['slope'].append(round(slope,4))
        euler_data['yn'].append(round(yn,4))

        print('------------------------------')
        y0 = yn
        x0 = x0+h
    
    print('\nAt x=%.4f, y=%.4f' %(xn,yn))
    result = (xn,yn)
    return result
##########   Euler Method Ends Here   ############

########## Trapezoidal Integration Method  #########
def trapezoid_f(func,val):
    x = val
    return 1/eval(func)
## Actual trapezoidal method logic
def trapezoidal(func,x0,xn,n):
    # calculating step size
    h = (xn - x0) / n
    
    # Finding sum 
    integration = trapezoid_f(func,x0) + trapezoid_f(func,xn)
    
    for i in range(1,n):
        k = x0 + i*h
        integration = integration + 2 * trapezoid_f(func,k)
    
    # Finding final integration value
    integration = integration * h/2
    
    return integration

#############  Simpson's 1/3 Rule Method starts here   #######
# Define function to integrate
def simpson13_f(func,value):
    x = value
    return 1/eval(func)
# Implementing Simpson's 1/3 
def simpson13(func,x0,xn,n):
    # calculating step size
    h = (xn - x0) / n
    
    # Finding sum 
    integration = simpson13_f(func, x0) + simpson13_f(func, xn)
    
    for i in range(1,n):
        k = x0 + i*h
        
        if i%2 == 0:
            integration = integration + 2 * simpson13_f(func, k)
        else:
            integration = integration + 4 * simpson13_f(func, k)
    
    # Finding final integration value
    integration = integration * h/3
    
    return integration

###########  SIMPSON's 3/8 RULE STARTSN HERE   #########
# Define function to integrate
def simpson38_f(func,value):
    x = value
    return 1/eval(func)

# Implementing Simpson's 3/8
def simpson38(func,x0,xn,n):
    # calculating step size
    h = (xn - x0) / n
    
    # Finding sum 
    integration = simpson38_f(func, x0) + simpson38_f(func, xn)
    
    for i in range(1,n):
        k = x0 + i*h
        
        if i%2 == 0:
            integration = integration + 2 * simpson38_f(func, k)
        else:
            integration = integration + 3 * simpson38_f(func, k)
    
    # Finding final integration value
    integration = integration * 3 * h / 8
    
    return integration

############  Runge Kutta 4 Method   ##################
rk4_data = {'iter':[] , 'x0':[] , 'y0':[] , 'yn':[]}
# function to be solved
def rk4_f(func,val_x,val_y):
    x = val_x
    y = val_y
    return eval(func)
# Actual logic and function for RK-4
def rk4(func,x0,y0,xn,n):
    
    # Calculating step size
    h = (xn-x0)/n
    
    print('\n--------SOLUTION--------')
    print('-------------------------')    
    print('x0\ty0\tyn')
    print('-------------------------')
    for i in range(n):
        k1 = h * (rk4_f(func, x0, y0))
        k2 = h * (rk4_f(func, (x0+h/2), (y0+k1/2)))
        k3 = h * (rk4_f(func, (x0+h/2), (y0+k2/2)))
        k4 = h * (rk4_f(func, (x0+h), (y0+k3)))
        k = (k1+2*k2+2*k3+k4)/6
        yn = y0 + k
        print('%.4f\t%.4f\t%.4f'% (x0,y0,yn) )
        print('-------------------------')
        rk4_data['iter'].append(i+1)
        rk4_data['x0'].append(round(x0,4))
        rk4_data['y0'].append(round(y0,4))
        rk4_data['yn'].append(round(yn,4))
        y0 = yn
        x0 = x0+h
    
    print('\nAt x=%.4f, y=%.4f' %(xn,yn))
    result = (xn,yn)
    return result

###########  Weddle's Rule Starts Here   ###########
def weddle_f(func,value):
    num = 1;
    x = value
    denom = float(eval(func));
 
    return num / denom;

# Actual weddle logic function
def weddle(func, a, b):
     
    # Find step size h
    h = (b - a) / 6;
     
    # To store the final sum
    sum = 0;
     
    # Find sum using Weedle's Formula
    sum = sum + (((3 * h) / 10) * (weddle_f(func, a)
            + weddle_f(func, a + 2 * h)
            + 5 * weddle_f(func, a + h)
            + 6 * weddle_f(func, a + 3 * h)
            + weddle_f(func, a + 4 * h)
            + 5 * weddle_f(func, a + 5 * h)
            + weddle_f(func, a + 6 * h)));
 
    # Return the final sum
    return sum;

###########   Boole's Rule Starts Here   ##############
# for the given value of x
def boole_f(func, value):
    x = value
    return (1 / eval(func))

# Actual Boole's Rule Logic
def BooleRule(func,a, b):
     
    # Number of intervals
    n = 4
 
    # Computing the step size
    h = ((b - a) / n)
    sum = 0
 
    # Substituing a = 0, b = 4 and h = 1
    bl = (7 * boole_f(func,a) + 32 * boole_f(func,a + h) + 12 *
        boole_f(func,a + 2 * h)+32 * boole_f(func,a + 3 * h)+7 *
        boole_f(func,a + 4 * h))* 2 * h / 45
 
    sum = sum + bl
    return sum

############   RK2 Method   ##############
rk2_data = {'iter':[] , 'x0':[] , 'y0':[] , 'slope':[], 'yn':[]}
# function to be solved
def rk2_f(expr,x,y):
    val=x
    val=y
    return eval(expr)

def rk2(func,x0,y0,xn,h):
    
    # Calculating step size
    n = int((xn-x0)/h)
    y = y0;
    
    print('\n-----------SOLUTION-----------')
    print('------------------------------')    
    print('x0\ty0\tslope\tyn')
    print('------------------------------')
    for i in range(1, n + 1):
        k1 = h * rk2_f(func,x0, y);
        k2 = h * rk2_f(func,x0 + 0.5 * h, y + 0.5 * k1);
        y = y + (1.0 / 6.0) * (k1 + 2 * k2);
        x0 = x0 + h;
        slope = (1.0/2.0)*(k1 + k2)
        yn = y + h * slope
        print('%.4f\t%.4f\t%0.4f\t%.4f'% (x0,y,slope,yn) )
        print('------------------------------')        
        rk2_data['iter'].append(i)
        rk2_data['x0'].append(round(x0,4))
        rk2_data['y0'].append(round(y0,4))
        rk2_data['slope'].append(round(slope,4))
        rk2_data['yn'].append(round(yn,4))
        y = yn
        x0 = x0+h
    
    print('\nAt x=%.4f, y=%.4f' %(xn,yn))
    res = [xn , yn]
    return res
################################################################################################
@app.route("/", methods=['GET' , 'POST'])
def index():
    bisection_data = {'iter':[],'a':[],'b':[],'c':[],'f(c)':[]}
    if request.method == 'POST':
        input_data = request.form.to_dict()
        #########   BISECTION METHOD   #########

        if(input_data['method'] == 'bisection'):
            #bisection_data = input_data
            function = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            eps = float(input_data['eps'])  #10E-6  # Acceptable Error (termination criteria)
            a = float(input_data['value1'])    #-2        # Guess Value for the lower bound on the root
            b = float(input_data['value2'])  #3        # Guess Value for the upper bound on the root

            if f(a,function)*f(b,function)>0:
                error = 'The given guesses do not bracket the root.'
                return render_template('index.html' , error=error)
            for i in range (maxiter):
                # Calculate the value of the root at the ith step
                c = (a+b)/2

                bisection_data['iter'].append(i+1)
                bisection_data['a'].append(round(a,4))
                bisection_data['b'].append(round(b,4))
                bisection_data['c'].append(round(c,4))
                bisection_data['f(c)'].append(round(f(c,function),4))
                #print(str(i+1)+'\t\t% 10.6f\t% 10.6f\t% 10.6f\t% 10.6f\t' %(a, b, c, f(c,function)))

                # Check if the root has been found with acceptable error or not?
                if np.abs(f(c,function))<eps:
                    root =c
                    #exit()
                # Check whether the root lies between a and c
                if f(a,function)*f(c,function)<0:
                    # Change the upper bound
                    b = c
                    root = b
                else: # The root lies between c and b
                    # Change the lower bound
                    a = c
                    root = a

                if i==maxiter-1:
                    print('\n\nMax iterations reached!')
                    print('Approximaiton to the Root after max iterations is : '+str(c))
            return render_template('index.html', bisection_input = input_data, bisection_data=bisection_data, bisection_root = root)
        
        ############  NEWTON RAPHSON METHOD   ######

        if(input_data['method'] == 'newtonRaphson'):
            expr = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            tolerable_err = float(input_data['tolerable_err'])
            x0 = float(input_data['guess'])

            #print(function,maxiter,tolerable_err,x0)
            #print("Expression : {} ".format(expr))

            deriv = derivative(expr)
            
            #print(newton_raphs_f(3,derivative(expr)))
            #print("Derivative of expression with respect to x : {}".format(expr_diff))
            #print("Value of the derivative : {} ".format(expr_diff.doit()))

            root = newtonRaphson(expr,x0,tolerable_err,maxiter)
            return render_template('index.html', newton_raphs_input = input_data, newton_raphs_data=newton_raphs_data, newton_raphs_root = root, derivative=deriv)

        ############   SECANT METHOD   ############
        if(input_data['method'] == 'secant'):
            expr = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            tolerable_err = float(input_data['tolerable_err'])
            x0 = float(input_data['guess1'])
            x1 = float(input_data['guess2'])

            root = secant(x0,x1,tolerable_err,maxiter,expr)
            return render_template('index.html/#secant', secant_input = input_data, secant_data=secant_data, secant_root = root)

        ############   REGULA FALSI METHOD   ############
        if(input_data['method'] == 'falsi'):
            expr = input_data['function']
            tolerable_err = float(input_data['tolerable_err'])
            x0 = float(input_data['guess1'])
            x1 = float(input_data['guess2'])
            root = falsi(x0,x1,tolerable_err,expr)
            print(falsi_data)
            return render_template('index.html', falsi_input =input_data, falsi_data=falsi_data, falsi_root = root)
        
        ############  SIMPLE ITERATIVE METHOD   ######

        if(input_data['method'] == 'fixedpoint'):
            expr = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            tolerable_err = float(input_data['tolerable_err'])
            x0 = float(input_data['guess'])

            root = fixedpoint(expr,x0,tolerable_err,maxiter)
            return render_template('index.html', fixedpoint_input = input_data, fixedpoint_data=fixedpoint_data, fixedpoint_root = root)
        ########    EULER METHOD STARTS   ##########
        if(input_data['method'] == 'euler'):
            func = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            xn = int(input_data['xn']) #calc point
            x0 = float(input_data['value1'])
            y0 = float(input_data['value2'])

            # Euler method call
            result = euler(func,x0,y0,xn,maxiter,input_data)
            xn = result[0]
            yn = result[1]
            return render_template('index.html',euler_data=euler_data, euler_input=input_data,xn=xn, yn=yn)
        ############   EULER METHOD ENDS    ##########
        
        ############  Trapexoidal Integ. Method  ##########
        if(input_data['method'] == 'trapezoidal'):            
            expr = input_data['function']
            x0 = float(input_data['x0']) # Upper Limit
            xn = float(input_data['xn']) # Lower Limit
            n = int(input_data['n']) # no. of sub intervals
            result = trapezoidal(expr, x0, xn, n)
            print (result)
            return render_template('index.html',trapezoidal_input=input_data, trapezoidal_result = result)

        ############   LANGRANGE METHOD STARTS   ###########
        if(input_data['method'] == 'langrange'):            
            n = int(input_data['datapoints'])
            xp = float(input_data['xp'])
            datapoints = input_data
            datapoints.pop('datapoints')
            datapoints.pop('method')
            datapoints.pop('xp')
            
            # Making numpy array of n & n x n size and initializing 
            # to zero for storing x and y value along with differences of y
            x = np.zeros((n))
            y = np.zeros((n))

            # separating values of x
            i = 0
            xi = 0
            yi = 0
            # print(datapoints['x0'])
            for key, value in datapoints.items():
                if((i % 2) == 0 ):
                    x[xi] = int(datapoints[key])
                    xi = xi + 1
                else:
                    y[yi] = int(datapoints[key])
                    yi = yi + 1
                i = i+1

            # Set interpolated value initially to zero
            yp = 0

            # Implementing Lagrange Interpolation
            for i in range(n):
                
                p = 1
                
                for j in range(n):
                    if i != j:
                        p = p * (xp - x[j])/(x[i] - x[j])
                
                yp = yp + p * y[i]    

            # Displaying output
            print('Interpolated value at %.3f is %.3f.' % (xp, yp))
            return render_template('index.html',langrange_input=input_data, xp=xp, yp=yp)
        ###########  LANGRANGE METHOD ENDS HERE   ########

        ###########  SIMPSON's 1/3 RULE STARTSN HERE   #########
        if(input_data['method'] == 'simpson13'):            
            expr = input_data['function']
            x0 = float(input_data['x0']) # Upper Limit
            xn = float(input_data['xn']) # Lower Limit
            n = int(input_data['n']) # no. of sub intervals
            result = simpson13(expr, x0, xn, n)
            print (result)
            return render_template('index.html',simpson13_input=input_data, simpson13_result = result)

        ###########  SIMPSON's 3/8 RULE STARTSN HERE   #########
        if(input_data['method'] == 'simpson38'):            
            expr = input_data['function']
            x0 = float(input_data['x0']) # Upper Limit
            xn = float(input_data['xn']) # Lower Limit
            n = int(input_data['n']) # no. of sub intervals
            result = simpson38(expr, x0, xn, n)
            print (result)
            return render_template('index.html',simpson38_input=input_data, simpson38_result = result)

        ########    RUNGE KUTTA 4 METHOD STARTS   ##########
        if(input_data['method'] == 'rk4'):
            func = input_data['function']
            maxiter = int(input_data['iter']) # Max. number of iterations
            xn = int(input_data['xn']) #calc point
            x0 = float(input_data['value1'])
            y0 = float(input_data['value2'])

            # Euler method call
            result = rk4(func,x0,y0,xn,maxiter)
            xn = result[0]
            yn = result[1]
            return render_template('index.html',rk4_data=rk4_data, rk4_input=input_data,xn=xn, yn=yn)
        
        ##########   WEDDLE's METHOD STARTS HERE   ############
        if(input_data['method'] == 'weddle'):
            func = input_data['function']
            a = float(input_data['a'])  # upper limit
            b = float(input_data['b'])  # lower limit

            result = round(weddle(func,a,b),4)
            #print(result)            
            return render_template('index.html',weddle_output=result , a=a, b=b, function=func)
        
        ##########   BOOLE'S METHOD STARTS HERE   ############
        if(input_data['method'] == 'boole'):
            func = input_data['function']
            a = float(input_data['a'])  # upper limit
            b = float(input_data['b'])  # lower limit

            result = round(BooleRule(func,a,b),4)
            #print(result)            
            return render_template('index.html',boole_output=result , a=a, b=b, function=func)
            
        ##############  RUNGE KUTTA 2 METHOD STARTS   ################
        if(input_data['method'] == 'rk2'):
            func = input_data['function']
            xn = int(input_data['xn']) #calc point
            x0 = float(input_data['value1'])
            y0 = float(input_data['value2'])
            h = float(input_data['h'])

            # Euler method call
            result = rk2(func,x0,y0,xn,h)
            xn = result[0]
            yn = result[1]
            return render_template('index.html',rk2_data=rk2_data, rk2_input=input_data,xn=xn, yn=yn)


    return render_template('index.html')

# Driver Code
if __name__ == "__main__":
    app.run(debug=True)