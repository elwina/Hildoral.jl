using Hildoral

A = [1 9
    1.1 10
    1.5 9.5
    1.4 8.5
    10.2 2
    9 1
    12 1.2
    12.5 2.2
    1.5 9.8
]
B = [0; 0; 0; 0; 1; 1; 1; 1; 1]
C = [2 5
    9 3]
ω, middle, jud, predre, err = Fisher(A, B)
Fisher2DPlot(A, B, ω, C)
