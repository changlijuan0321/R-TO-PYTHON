#-*- coding:utf-8 -*-
import json
import codecs

path = '47808-bacnet-device_id-full_ipv4-20170623T020002-zgrab-results.json'

with open(path, 'r', encoding='utf-8') as f:
    with open('ip-location.txt', 'w', encoding='utf-8') as tmp:
        count = 0
        s = f.readline()
        while s:
            try:
                new_dict = json.loads(s)
                count += 1
                # 读的这一行，如果 location 有值就把 location 对应的键值和对应的
                # ip 一起存入一个文件
                if ("location" in new_dict["data"]["bacnet"].keys() 
                and new_dict["data"]["bacnet"]["location"] != "unknown"):
                    tmp.write((new_dict["ip"] + new_dict["data"]["bacnet"]["location"]))
                s = f.readline()
            except Exception as e:
                print(e)
                continue
        print(count)

        
        

                




