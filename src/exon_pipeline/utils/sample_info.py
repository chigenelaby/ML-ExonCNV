#-*-coding:utf-8-*-
"""
返回样本信息类
"""

__version__ = '0.1.0'

import os
import configparser


class WkCodeBase:
    """文库号元类
    提供公共接口
        wkcode  文库号
        analyse_path    分析目录
        path      分析目录
    
    >>> wkcode = WkCodeBase('AAV188', 'AAV188_NTb01F_NT01T_b53591')
    >>> wkcode 
    AAV188
    >>> wkcode == 'AAV188'
    True
    """
    def __init__(self, wkcode, analyse_path=None, config_file=None):
        self.wkcode = wkcode
        self.analyse_path = analyse_path
        self.config = None
        if config_file:
            self.read_config(config_file)
    
    def read_config(self, config_file):
        """读取config文件"""
        config = configparser.ConfigParser()
        config.read(config_file)
        self.set_config(config)
    
    def set_config(self, config):
        """配置config"""
        self.config = config
    
    def get_filename(self, file):
        """获取配置文件中文件标识名"""
        
        if file in self.config.options('base'):
            return self.config['base'][file]
        if file in self.config.options('filename'):
            filename = self.config['filename'][file]
            return os.path.join(self.path, filename)
    
    @property
    def path(self):
        return self.analyse_path
    
    def __str__(self):
        return self.wkcode
    
    def __repr__(self):
        return self.wkcode
    
    def __eq__(self, other):
        if hasattr(other, 'wkcode'):
            return self.wkcode == other.wkcode
        if isinstance(other, str):
            return self.wkcode == other
        return False
    
    def __hash__(self):
        """重写hash"""
        return hash(self.wkcode)
    
    def forwark(self):
        """初始化"""

def get_wkcode(wkcode, analyse_path, cls=WkCodeBase, config_file=None, config=None):
    """获取文库号类"""
    wk = cls(wkcode, analyse_path, config_file)
    if not config_file and config:
        wk.set_config(config)
    return wk
