
from optparse import OptionParser
import solver
import plotter
from sys import argv, exit
 
if __name__ == '__main__':
        parser = OptionParser()
        parser.add_option('-n', dest='n', default=10)
        parser.add_option('-o', '--output', dest='output_file', default='solution.txt')
        
        (options, args) = parser.parse_args()
        options.n = int(options.n)

        (A, b, x) = solver.to_linear_system(options.n) 
        print(A.size, b.size, x.size)
        y = solver.solve_system(A, b)
        print(y)
        print(x)

        plotter.plot(y, x)

        print('Start solving the model...')
        print(b)
        print(A)
        print('Done!')
        
        f = open(options.output_file, 'w')
        f.close()
        
        
