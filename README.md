# JPEG_LS_ARGB:一种ARGB数据无损压缩/解压算法
*Contributor: Ausont Lu, Little Chestnut*


## 1 项目介绍

**项目背景**

本项目脱胎于第七届全国大学生集成电路创新竞赛的景嘉微杯赛题，杯赛题目为杯赛题目：**一种ARGB数据无损压缩/解压单元**，目前开源的无损压缩算法仍然不充分、代码和算法原理对应的可读性仍然不清晰，对于初学者较不友好，因此将我们在竞赛中所做出的算法开源。

**算法原理**

JPEG-LS是一种经典的无损压缩算法，但该算法最初被应用于单通道的灰度图像。我们基于国际电联T.87标准的JPEG-LS算法，将其应用于ARGB的四通道图像处理，方法是将四个通道并行的分别通过JPEG-LS进行压缩/解压。

**代码贡献**

本项目基于景嘉微公司提供的压缩/解压程序接口开发，程序外层的文件流处理、MakeFile由景嘉微公司开发，核心JPEG-LS算法参考了davidebardone的Github项目，但该项目未完成游程编解码这一核心功能，并有一定的bug。我们在该项目的基础上完成了游程编解码的功能，修复了bug并拓展到了ARGB四个通道，并将代码模块按T.87标准中的模块严格对应。


## 2 编译与运行
### 2.1 开发运行环境
| 项目     | 要求|
| -------- | ------------------------------------------------------------ |
| 操作系统 | ubuntu 18或以上 |
| 开发语言 | C/C++，如使用C语言标准库以外的其他软件，则该软件源代码文件必须包含在程序中 |
| 其他要求 | 使用Makefile编译 |

### 2.2 程序接口与使用方法

**可执行文件名**  
`fblcd.out`（其中：`fb` 含义为 framebuffer；`lcd`含义为 Lossless Compression and Decompression）

**命令行参数**

命令行格式为： 

```
fblcd.out [--version] [-{en,de} infile outfile]
```

参数含义：

| 参数名称  | 含义 |
| --------- | ------------------------------------------------------------ |
| --version | 显示程序版本号等信息 |
| -en       | 执行ARGB数据压缩。<br />*infile* 指定输入文件的路径，输入文件为BMP格式的图像文件；*outfile* 指定输出文件的路径。<br />*infile*、*outfile* 的路径可以是相对路径，也可以是绝对路径 |
| -de       | 执行数据解压缩，得到ARGB数据。<br />*infile*、*outfile* 的含义与 -en 的相同。其中 *outfile* 所指定输出文件为BMP格式的图像文件 |

当程序运行时没有指定任何参数，则显示命令行格式，并返回-1。

**程序返回值含义**

程序主函数`main()`的返回值用于表示执行成功与否，定义如下：

| 返回值 | 含义               |
| ------ | ------------------ |
| 0      | 执行成功           |
| -1     | 参数不足           |
| -2     | 无效参数           |
| -3     | 输入文件打开失败   |
| -4     | 输出文件创建失败   |
| -5     | 输入文件格式不正确 |
| 其他   | （请自行定义）     |

### 2.3 输入、输出数据格式

**BMP格式图像文件**

数据压缩的输入文件与数据解压缩的输出文件为BMP格式的图像文件，每个像素的位数为32位，图像的宽度和高度都必须是8的倍数。
【注意】样本数据文件中是32位的，每一个颜色通道中都有有效数据，请不要做文件格式转化或编辑，否则Alpha数据可能丢失。

**压缩格式数据文件**

数据压缩的输出文件与数据解压缩的输入文件为压缩格式文件，其文件格式如下：

| 数据块                | 偏移量 | 字节数        | 内容                                                         |
| --------------------- | ------ | ------------- | ------------------------------------------------------------ |
| 文件标识              | 0000H  | 4             | 固定为4个字符：`JLCD`                                        |
| Width                 | 0004H  | 4             | 图像宽度                                                     |
| Height                | 0008H  | 4             | 图像高度                                                     |
| TileWidth             | 000CH  | 4             | Tile宽度，应当是8或16                                        |
| TileHeight            | 0010H  | 4             | Tile高度，应当是8或16                                        |
| TileCount             | 0014H  | 4             | Tile总数，其值应当为 Width * Height / (8 * 8)                |
| DataInfo              | 0018H  | 8 * TileCount | 从Tile[0]开始的所有Tile的信息，每个Tile有2个32位整数：<br />（1）起始偏移地址（TilePos），4 Bytes。<br />偏移地址含义为：对于Tile[i]，其值为`TileData[i]在文件中的位置 - TileData[0]在文件中的位置`<br />（2）该Tile的长度（TileLen），4 Bytes，以字节为单位 |
| TileData[0]           | x      | x             | Tile[0]数据                                                  |
| TileData[1]           | x      | x             | Tile[1]数据                                                  |
| ……                    | x      | x             |                                                              |
| TileData[TileCount-1] | x      | x             | 最后一个Tile的数据                                           |

【说明】以上数据中，整数数据使用小端字节序，即`0x12345678`在内存中的数据为 `0x78 0x56 0x34 0x12`。

## 3 JPEG-LS算法核心函数部分
算法核心代码位于rgbTileProc.h及rgbTileProc.cpp文件中，其余文件处理文件流及输入输出。其中，负责编码的核心函数为argb2tile，解码的核心函数为tile2argb，两个函数输入输出变量的具体含义可见代码文件的注释。

## 4 参考资料

**[1] ITU-T T.87** : Information technology – Lossless and near-lossless compression of continuous-tone still images – Baseline :  https://www.itu.int/rec/T-REC-T.87/en

**[2] davidebardone的Github项目 jpeg-ls**  : https://github.com/davidebardone/jpeg-ls
