## IMDB : Infectious diseases pathogens' Mutation profile Database

IMDB是一个交互式在线Web应用平台，用于分析呼吸道传染病样本中的突变，并提供基于突变谱的生物应用程序，如自动qPCR检测设计。

### 部署

#### 0. 环境要求：

python >= 3.6<br>
uwsgi

django<br>
djangorestframework<br>
markdown<br>
django-filter<br>
pybind11<br>
pandas<br>
biopython<br>
django-cors-headers<br>
mysqlclient<br>
pymysql

#### 1. 下载

git clone https://github.com/catsingchannel/IMDB.git

#### 2. 编译C++组件

##### pysnp：

进入pysnp文件夹，然后运行以下命令：

```
mkdir build
cd build
cmake ..
make
```

复制build文件夹中的pysnpprocess.cpython-3XX-x86_64-linux-gnu.so文件到/项目路径/app目录下

##### primerDesign：

进入primerDesign文件夹，然后运行以下命令：

```
mkdir build
cd build
cmake ..
make
```

复制build文件夹中的primer_design.cpython-3XX-x86_64-linux-gnu.so文件到/项目路径/app目录下

#### 3. 初始化数据库

##### 若使用默认数据库则跳过此步骤

进入项目文件夹内运行：
```
python3 manage.py migrate
```

#### 4. uwsgi启动服务

新建一个uwsgi配置文件（后缀为ini）

参数配置参考：
```
socket = 127.0.0.1:[端口号]
chdir = /路径/项目文件夹
module = GVMDF.wsgi:application
master = true
workers = 2
reload-mercy = 60
vacuum = true
http-timeout = 7200
harakiri = 7200
pidfile= /路径/pid文件名
daemonize= /路径/日志文件名
```

然后启动uwsgi服务：

```
uwsgi --ini /路径/配置文件名.ini
```

#### 5. nginx代理设置

api接口参考设置：
```
    location /api/ {
	    rewrite ^/api/(.*)$ /$1 break;
	    proxy_pass http://127.0.0.1:[端口号];
	    proxy_set_header Host $host;
	    proxy_set_header X-Real-IP $remote_addr;
	    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
	    proxy_set_header X-Forwarded-Proto $scheme;
	    proxy_read_timeout 7200s;
	    proxy_connect_timeout 7200s;
    }

```

静态文件接口参考设置：
```
location /static/ {
	    alias /路径/项目文件夹/static/;
    }

```


#### 6. 访问

前端页面为/app/templates/homepage.html<br>
JavaScript程序为/app/JavaScript.js

使用时修改这两个文件内的域名为对应api域名，访问时发送homepage.html文件即可

### 上传数据

#### 0. 环境要求与文件格式

##### 环境要求：

MUMmer 4.0

##### 文件格式：

上传数据需要两个文件，一个fasta格式保存样本序列，一个csv格式保存样本信息<br>
文件格式参见./app/upload文件夹下的两个example文件<br>
若对应物种基因组不是分段结构，则csv文件需要去除“segment”列，或者将该列值全部赋为0。

#### 1. 数据上传

##### 0. 将fasta和csv两个文件放到 /项目路径/app/upload/ 文件夹下

##### 1. 启动django shell：

在项目文件夹下启动django shell:
```
python3 manage.py shell
```

##### 2. 上传数据：

```
from app.process import data_entry

s = "你的物种名称"
#以example中物种为例的话是 'human Influenza A Virus'
#所有物种名参见网页中Species一栏中的选项
#或数据库中查询app_speceis表/ORM 中查询Species模型

data_entry(s)
```

数据上传时间取决于上传样本量的大小，大量样本的上传时间会较为缓慢且会占用大量内存，建议分割上传
