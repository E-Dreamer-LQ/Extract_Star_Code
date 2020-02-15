# Extract_Star_Code

提星流程

    1. 首先找到fits文件，通过加载文件路径完成。
    2. 找到fits文件后通过open()和for循环迭代完成对fits数据的读取。读出来的数据类型为numpy.ndarray。
    3. 算背景包括将图像分解为一个局部区域，计算其像素平均值，中值，标准差，设置一个阈值，将背景从图像中去除。
    4. 接着前面形成的局部区域和设置的阈值，通过扫描像素区域与阈值进行比较，做二值化处理，也就是形成一个只有0和1的矩阵。
    5. 得到的二值矩阵，通过bwlabel算法计算连通域。
    6. 依据噪声和星的差别，噪声是独立存在的即单像素存在，星是一块亮的像素区域即是连通的，所以需要一个去噪的处理，程序中是将存在的噪声位置坐标全部写入到txt文件。

