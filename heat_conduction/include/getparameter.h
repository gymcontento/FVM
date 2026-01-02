#ifndef GETPARAMETER_H
#define GETPARAMETER_H
#include <string>
#include <unordered_map>
#include <variant>
#include <fstream>
#include <iostream>
#include <sstream>

using ParameterValue = std::variant<std::string, int, float, double>;

class Getparameter
{
public:
    virtual ~Getparameter() = default;
    virtual bool exportParameters(std::string& filename) = 0;
};

class GeoJsonParameter : public Getparameter
{
private:
    std::unordered_map<std::string, ParameterValue> parameters;
public:
    bool exportParameters(std::string& filename)
    {
        std::ifstream _ifs(filename);
        if(!_ifs.is_open())
        {
            std::cerr << "Open parameter File failed!" << std::endl;
            return false;
        }

        // 读取整个文件内容
        std::string content((std::istreambuf_iterator<char>(_ifs)),
                        std::istreambuf_iterator<char>());

        // 移除所有空白字符
        content.erase(std::remove_if(content.begin(), content.end(), ::isspace), content.end());
        
        // 移除开头和结尾的花括号
        if(content.front() == '{') content.erase(0, 1);
        if(content.back() == '}') content.pop_back();
        
        // 按逗号分割键值对
        std::stringstream ss(content);
        std::string pair;

        while(std::getline(ss, pair, ',')) {
            // 找到冒号位置
            size_t colon_pos = pair.find(':');
            if(colon_pos == std::string::npos) continue;
            
            // 提取键
            std::string key = pair.substr(0, colon_pos);
            // 提取值
            std::string value = pair.substr(colon_pos + 1);
            
            // 移除键的引号
            if(key.front() == '"' && key.back() == '"') {
                key = key.substr(1, key.length() - 2);
            }
            
            // 移除值的引号
            if(value.front() == '"' && value.back() == '"') {
                value = value.substr(1, value.length() - 2);
            }
            
            // 尝试判断值的类型并存储
            try {
                // 尝试转换为double（可以处理整数和浮点数）
                double double_val = std::stod(value);
                
                // 检查是否为整数
                if(value.find('.') == std::string::npos) {
                    parameters[key] = static_cast<int>(double_val);
                } else {
                    parameters[key] = double_val;
                }
            } catch(...) {
                // 如果转换失败，作为字符串存储
                parameters[key] = value;
            }
        }
        return true;
    }
};

#endif

