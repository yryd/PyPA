# database.py
# 模拟反应条件与结果数据库

import logging
import sqlite3
import csv
import os
import math
from src.generator import generate_reaction_smile
from src.molecular import MolecularModule

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# 数据库模块类
class DatabaseModule:
    def __init__(self, db_path='data/database/reaction.db'):
        """初始化数据库模块，创建数据库连接和游标
        Args:
            db_path: 数据库文件的路径（包含文件名）
        """
        self.path = db_path
        self.r1_list = []
        self.r2_list = []
        self.sol_list = []
        self.rat_list = []
        self.connection = sqlite3.connect(db_path)
        self.cursor = self.connection.cursor()

    def data_init(self, base_num : int = 250, properties_list : list = [
        'group_smiles', 'rings', 'count_group', 
        'group_min_distances', 'group_max_distances', 'ar_sp3_balance', 
        'aromatic_rings', 'molfinger'
    ]):
        """初始化数据库并填充数据

        该方法执行一系列操作以初始化数据库，包括创建数据表、生成反应产物并存储、提取分子的结构特征并记录到数据库中，
        最后将数据库导出为CSV文件。

        Args:
            base_num (int): 反应规模，体系中分子总数，默认为250。
            properties_list (list): 要提取并存储的分子属性名称列表，默认为包含所有分子特征的列表。
                                    这些属性将被提取并添加到数据库中。
        """
        try:
            self.create_table()
            logging.info("创建反应列表完成")
            # 根据反应产生产物 SMILES 字符串，填充数据库
            logging.info("等待生成产物......")
            self.generate_and_store_product()
            logging.info("生成产物 SMILES 字符串、构建反应索引完成")
            # 对分子提取某些结构特征记入数据库
            logging.info("等待填充反应物分子特征属性信息...")
            self.fill_data_properties(properties_list)
            logging.info("填充分子特征属性完成")
            # 提取反应信息：分子模拟分子个数
            self.add_molecular_num(base_num)
            logging.info("计算体系分子数量完成")
            # 数据库导出为CSV
            self.export_to_csv()
            logging.info("成功初始化数据库，可查看csv文件校对数据库！")
        except Exception as e:
            logging.error(f"构建lt文件出错: {e}")
            raise

    def create_table(self, smile_path='data/smiles/'):
        """创建初始化反应物条件数据表并填充数据
        Args:
            smile_path: SMILE文件的路径（仅路径）
        """
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS reactions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                reactant1_smiles TEXT,
                reactant2_smiles TEXT,
                solvent_smiles TEXT,
                reactant1_ratio INTEGER,
                reactant2_ratio INTEGER,
                solvent_ratio INTEGER,
                reactant1_key INTEGER,
                reactant2_key INTEGER,
                solvent_key INTEGER
            )
        ''')
        self.connection.commit()
        
        all_list = self.get_reactant_list(smile_path)
        # 添加反应物基础 SMILE 信息
        self.fill_data(*all_list)

    def get_smiles_file(self, file_path):
        """读取SMILES文件并返回每行的列表
        Args:
            file_path: SMILES文件的路径
        Returns:
            返回包含SMILES字符串的列表
        """
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        # 去掉每行末尾的换行符，并忽略空行
        lines = [line.strip() for line in lines if line.strip()]
        return lines

    def get_ratios(self):
        """获得反应物、溶剂的比例，两种反应物比例不为0。
        Returns:
            返回包含可能反应物比例的元组列表
        """
        results = []
        for a in range(1, 5):
            for b in range(1, 5):
                c = 5 - a - b
                if c >= 0:
                    results.append((a, b, c))
        return results

    def get_reactant_list(self, smile_file_path):
        """获取反应物和溶剂的SMILES列表
        Args:
            smile_file_path: 反应物与溶剂SMILES文件的路径
        Returns:
            返回包含反应物1、反应物2、溶剂和比例的元组
        """
        reactant1_list = self.get_smiles_file(smile_file_path + 'smile1.list')
        reactant2_list = self.get_smiles_file(smile_file_path + 'smile2.list')
        solvent_list = self.get_smiles_file(smile_file_path + 'smile3.list')
        ratios_list = self.get_ratios()
        return (reactant1_list, reactant2_list, solvent_list, ratios_list)

    def fill_data(self, reactant1_list, reactant2_list, solvent_list, ratios):
        """通过遍历列表初步填充数据
        Args:
            reactant1_list: 反应物1的SMILES字符串列表
            reactant2_list: 反应物2的SMILES字符串列表
            solvent_list: 溶剂的SMILES字符串列表
            ratios: 包含反应物和溶剂比例的元组列表
        """
        # 遍历反应物和溶剂，插入数据
        for idx1, r1 in enumerate(reactant1_list):
            for idx2, r2 in enumerate(reactant2_list):
                for idx_solvent, solvent in enumerate(solvent_list):
                    for ratio in ratios:
                        r1_ratio, r2_ratio, solvent_ratio = ratio

                        # 插入数据与分子键值
                        self.cursor.execute('''
                            INSERT INTO reactions (reactant1_smiles, reactant2_smiles, solvent_smiles, 
                                                reactant1_ratio, reactant2_ratio, solvent_ratio,
                                                reactant1_key, reactant2_key, solvent_key)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (r1, r2, solvent, r1_ratio, r2_ratio, solvent_ratio, idx1, idx2, idx_solvent))

        self.connection.commit()

    def read_data(self):
        """读取数据库中的所有数据
        Returns:
            包含数据库数据的元组列表
        """
        self.cursor.execute('SELECT * FROM reactions')
        result = self.cursor.fetchall()
        return result

    def add_column(self, column_name, column_type):
        """在表中添加新列
        Args:
            column_name: 新列的名称
            column_type: 新列的数据类型（如REAL, TEXT等）
        """
        self.cursor.execute(f'ALTER TABLE reactions ADD COLUMN {column_name} {column_type}')
        self.connection.commit()

    def check_and_add_columns(self, columns):
        """根据列名检查列是否存在，并添加缺失的列
        
        Args:
            columns (dict): 包含列名及其类型的字典，例如 {'product_smiles': 'TEXT', 'byproduct_smiles': 'TEXT'}
        """
        # 获取当前表的列名
        self.cursor.execute("PRAGMA table_info(reactions)")
        existing_columns = {column[1] for column in self.cursor.fetchall()}

        for column_name, column_type in columns.items():
            if column_name not in existing_columns:
                self.add_column(column_name, column_type)

    def add_data_column(self, column_name, input_list, column_type = 'REAL'):
        """根据列名添加新数据列，并根据输入列表添加新行，为所有原数据添加新特征。
        此功能可根据列表将数据量翻倍，如果输入长度为1的列表则为所有数据添加一个新特征，
        如果输入列表长度为2，则新特征的不同取值与原数据正交，数据量翻一倍，更长的列表以此类推
        Args:
            column_name: 新列的名称
            input_list: 新的数据列表
        """
        # 添加新列
        self.add_column(column_name, column_type)

        # 读取现有数据
        existing_data = self.read_data()

        # 获取列名列表（不包括 id 列）
        column_names = [column[0] for column in self.cursor.description if column[0] != 'id']

        # 删除原数据库中的所有内容
        self.cursor.execute('DELETE FROM reactions')
        self.connection.commit()

        # 重置主键
        self.cursor.execute('DELETE FROM sqlite_sequence WHERE name="reactions"')
        self.connection.commit()

        # 遍历现有数据并添加新行
        for row in existing_data:
            for value in input_list:
                new_row = list(row[1:-1]) + [value]  # 复制现有行（去掉 id 列）并添加新列的值
                placeholders = ', '.join(['?'] * len(new_row))  # 创建占位符
                self.cursor.execute(f'''
                    INSERT INTO reactions ({', '.join(column_names)})
                    VALUES ({placeholders})
                ''', new_row)
        
        self.connection.commit()

    def get_column_list(self, column_name):
        """根据列名获取指定列的所有值
        Args:
            column_name: 要查询的列名
        Returns:
            包含该列所有值的列表
        """
        self.cursor.execute(f'SELECT {column_name} FROM reactions')
        column_values = [row[0] for row in self.cursor.fetchall()]  # 获取所有列值
        return column_values

    def get_row_by_id(self, row_id):
        """根据id获取整行数据
        Args:
            row_id: 要查询的行的id
        Returns:
            包含该行数据的列表，如果没有找到则返回 None
        """
        self.cursor.execute('SELECT * FROM reactions WHERE id = ?', (row_id,))
        row = self.cursor.fetchone()  # 获取单行数据
        return list(row) if row else None  # 返回行数据的列表或 None

    def get_value_by_id(self, row_id, column_name):
        """根据id与列名获取单个值
        Args:
            row_id: 要查询的行的id
            column_name: 要获取的列名
        Returns:
            指定列的值，如果没有找到则返回 None
        """
        query = f'SELECT {column_name} FROM reactions WHERE id = ?'
        self.cursor.execute(query, (row_id,))
        row = self.cursor.fetchone()  # 获取单行数据

        return row[0] if row else None  # 返回指定列的值或 None

    def update_value(self, row_id, column_name, new_value):
        """根据id和列名更新指定列的值
        Args:
            row_id: 要更新的行的id
            column_name: 要更新的列名
            new_value: 新的值
        """
        self.cursor.execute(f'''
            UPDATE reactions
            SET {column_name} = ?
            WHERE id = ?
        ''', (new_value, row_id))
        self.connection.commit()

    def export_to_csv(self, csv_file_path = 'data/database/reaction.csv'):
        """将数据库导出为CSV文件
        Args:
            csv_file_path: 导出CSV文件的路径（包含文件名）
        """
        data = self.read_data()
        with open(csv_file_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            # 写入列名
            writer.writerow([column[0] for column in self.cursor.description])
            # 写入数据
            writer.writerows(data)

    def delete_database(self):
        """删除数据库文件
        
        Raises:
            FileNotFoundError: 如果数据库文件未找到
            OSError: 如果删除数据库时发生其他错误
        """
        try:
            os.remove(self.path)
            print(f"数据库 {self.path} 已删除。")
        except FileNotFoundError:
            raise FileNotFoundError(f"数据库文件 {self.path} 未找到。")
        except OSError as e:
            raise OSError(f"删除数据库时发生错误: {e}")

    def close(self):
        """关闭数据库连接"""
        self.connection.close()
    
    def generate_and_store_product(self):
        """为所有反应生成反应产物与反应索引并存储"""
        
        # 检查并添加缺失的列
        self.check_and_add_columns({
            'product_smiles': 'TEXT',
            'byproduct_smiles': 'TEXT',
            'reaction_index_dicts': 'TEXT'
        })
        self.cursor.execute('SELECT id, reactant1_smiles, reactant2_smiles FROM reactions')
        reactions = self.cursor.fetchall()

        # 创建一个缓存字典
        reaction_cache = {}

        for reaction in reactions:
            id, reactant1_smile, reactant2_smile = reaction
            # print(reaction)
            # # 生成反应产物
            # product_smiles_list, byproduct_smiles, atom_indices_dict_list = generate_reaction_smile(reactant1_smile, reactant2_smile)
            
            # 创建一个唯一的键用于缓存
            reaction_key = reactant1_smile + reactant2_smile
            
            # 检查缓存中是否已有结果
            if reaction_key not in reaction_cache:
                # 生成反应产物
                product_smiles_list, byproduct_smiles, atom_indices_dict_list = generate_reaction_smile(reactant1_smile, reactant2_smile)
                # 将结果存入缓存
                reaction_cache[reaction_key] = (product_smiles_list, byproduct_smiles, atom_indices_dict_list)
            else:
                # 从缓存中获取结果
                product_smiles_list, byproduct_smiles, atom_indices_dict_list = reaction_cache[reaction_key]
                
            # 将生成物和副产物存储到数据库，列表转换为字符串，使用 ; 分隔
            product_smiles = ';'.join(product_smiles_list)
            reaction_index_str = encode_nested_structure_v2(atom_indices_dict_list)
            
            self.cursor.execute('''
                UPDATE reactions
                SET product_smiles = ?, byproduct_smiles = ?, reaction_index_dicts = ?
                WHERE id = ?
            ''', (product_smiles, byproduct_smiles, reaction_index_str, id))

            self.connection.commit()

    def fill_data_properties(self, properties):
        """通过遍历列表初步填充数据
        Args:
            properties: 需要存储的分子属性名列表
        """
        # 获取列名
        self.cursor.execute("PRAGMA table_info(reactions)")
        existing_columns = {column[1] for column in self.cursor.fetchall()}
        # 执行查询以获取所有反应数据
        self.cursor.execute('SELECT id, reactant1_smiles, reactant2_smiles, solvent_smiles FROM reactions')
        # 获取所有行
        rows = self.cursor.fetchall()
        for row in rows:
            id, reactant1_smile, reactant2_smile, solvent_smile = row
            mol_r1 = MolecularModule(reactant1_smile, 'r1')
            mol_r1.cal_mol_prop()
            mol_r2 = MolecularModule(reactant2_smile, 'r2')
            mol_r2.cal_mol_prop()
            mol_solvent = MolecularModule(solvent_smile, 'sol')
            mol_solvent.cal_mol_prop()
            # 处理分子属性
            for prop in properties:
                # 先假定类型为 'REAL'，然后根据实际情况调整
                column_type = 'REAL'
                # 对属性添加前缀
                for mol, prefix in zip([mol_r1, mol_r2, mol_solvent], ['reactant1_', 'reactant2_', 'solvent_']):
                    if hasattr(mol, prop):
                        value = getattr(mol, prop)
                        # 处理不同类型的属性
                        if isinstance(value, list):
                            # 如果是列表，转换为字符串
                            value = encode_nested_structure(value)
                            column_type = 'TEXT'  # 列表类型
                        elif isinstance(value, str):
                            column_type = 'TEXT'  # 字符串类型
                        elif isinstance(value, int):
                            column_type = 'INTEGER'  # 整数类型
                        elif isinstance(value, float):
                            column_type = 'REAL'  # 浮点数类型


                        # 确保只为每个 prop 调用一次 add_column
                        column_name = f'{prefix}{prop}'
                        if column_name not in existing_columns:
                            self.add_column(column_name, column_type)
                            # 更新列名列表
                            self.cursor.execute("PRAGMA table_info(reactions)")
                            existing_columns = {column[1] for column in self.cursor.fetchall()}

                        self.update_value(id, column_name, value)
        return

    def add_molecular_num(self, base_num = 250):
        """创建新列，根据环数与比例判断模拟分子个数，再计算出反应物反应基团总数，并填充数据库
        该方法通过计算反应物和溶剂的比例、环数以及基数，来确定每种分子的数量和反应基团的数量，并将这些值更新到数据库中。

        Args:
            base_num (int): 总的基础分子数量，默认为250。用于计算每种分子的数量。
        """
        # 检查并添加新列，确保不重复
        self.check_and_add_columns({
            'r1_num': 'INTEGER',
            'r2_num': 'INTEGER',
            'sol_num': 'INTEGER',
            'r1_group_num': 'INTEGER',
            'r2_group_num': 'INTEGER'
        })
        # 获取所有反应数据
        self.cursor.execute('SELECT id, reactant1_ratio, reactant2_ratio, solvent_ratio, '
                            'reactant1_rings, reactant2_rings, '
                            'reactant1_count_group, reactant2_count_group FROM reactions')
        rows = self.cursor.fetchall()  # 获取所有行

        for row in rows:
            id, r1_ratio, r2_ratio, sol_ratio, r1_rings, r2_rings, r1_count_group, r2_count_group = row
            
            # 计算分子数量，向上取整
            r1_num = math.ceil(base_num * r1_ratio / (r1_ratio + r2_ratio + sol_ratio) / r1_rings)
            r2_num = math.ceil(base_num * r2_ratio / (r1_ratio + r2_ratio + sol_ratio) / r2_rings)
            sol_num = math.ceil(base_num * sol_ratio / (r1_ratio + r2_ratio + sol_ratio))
            
            # 计算反应基团数量
            r1_group_num = r1_num * r1_count_group
            r2_group_num = r2_num * r2_count_group
            
            # 更新数据库中的值
            self.cursor.execute('''
                UPDATE reactions
                SET r1_num = ?, r2_num = ?, sol_num = ?, 
                    r1_group_num = ?, r2_group_num = ?
                WHERE id = ?
            ''', (r1_num, r2_num, sol_num, r1_group_num, r2_group_num, id))
            
        # 提交更改
        self.connection.commit()
        return

    def add_product_fun_group(self):
        """获取产物分子官能团索引并填充到数据库中"""
        
        # 检查并添加新列
        self.check_and_add_columns({
            'product_mol_functional_group': 'TEXT'
        })

        # 获取所有 product_smiles
        self.cursor.execute('SELECT id, product_smiles FROM reactions')
        rows = self.cursor.fetchall()

        for row in rows:
            id, product_smiles = row
            
            # 将 product_smiles 分割成列表
            smiles_list = product_smiles.split(';')
            
            # 存储每个分子属性的列表
            property_values = []

            for smile in smiles_list:
                mol = MolecularModule(smile, 'p')
                property_value = getattr(mol, 'mol_functional_group')
                property_values.append(property_value)

            # 将属性值转换为字符串格式
            property_values_str = encode_nested_structure(property_values)

            # 更新数据库中的指定属性列
            self.cursor.execute(f'''
                UPDATE reactions
                SET 'product_mol_functional_group' = ?
                WHERE id = ?
            ''', (property_values_str, id))

        # 提交更改
        self.connection.commit()


def encode_nested_structure(nested):
    """
    将嵌套列表或元组编码为字符串。

    Args:
        nested (list or tuple): 嵌套列表或元组。

    Returns:
        str: 编码后的字符串。
    """
    if isinstance(nested, list):
        return '[' + ', '.join([encode_nested_structure(item) for item in nested]) + ']'
    elif isinstance(nested, tuple):
        return '(' + ', '.join(map(str, nested)) + ')'
    else:
        return str(nested)

def decode_nested_structure(encoded_str):
    """
    将编码后的字符串解码为嵌套列表或元组。

    Args:
        encoded_str (str): 编码后的字符串。

    Returns:
        list or tuple: 解码后的嵌套列表或元组。
    """
    # 去掉外层的方括号或圆括号
    if encoded_str.startswith('[') and encoded_str.endswith(']'):
        content = encoded_str[1:-1]
        return [decode_nested_structure(item.strip()) for item in split_nested(content)]
    elif encoded_str.startswith('(') and encoded_str.endswith(')'):
        content = encoded_str[1:-1]
        return tuple(decode_nested_structure(item.strip()) for item in split_nested(content))
    else:
        return eval(encoded_str)  # 直接使用 eval 处理基础数据类型

def split_nested(s):
    """
    将字符串按逗号分割，同时考虑嵌套结构。

    Args:
        s (str): 输入字符串。

    Returns:
        list: 分割后的字符串列表。
    """
    result = []
    current = []
    depth = 0

    for char in s:
        if char == '[' or char == '(':
            depth += 1
        elif char == ']' or char == ')':
            depth -= 1
        
        if char == ',' and depth == 0:
            result.append(''.join(current).strip())
            current = []
        else:
            current.append(char)

    if current:
        result.append(''.join(current).strip())

    return result

def encode_nested_structure_v2(nested):
    """
    将嵌套列表或字典编码为字符串。

    Args:
        nested (list or dict): 嵌套列表或字典。

    Returns:
        str: 编码后的字符串。
    """
    if isinstance(nested, list):
        return '[' + ', '.join([encode_nested_structure_v2(item) for item in nested]) + ']'
    elif isinstance(nested, dict):
        items = [f"{repr(key)}: {encode_nested_structure_v2(value)}" for key, value in nested.items()]
        return '{' + ', '.join(items) + '}'
    else:
        return repr(nested)


if __name__ == "__main__":
    db = DatabaseModule(db_path = '../data/database/reaction.db')
    db.create_table(smile_path='../data/smiles/')
    # db.add_data_column('new_line',[5555,321])
    print(db.get_row_by_id(5))
    # print(db.get_column_list('id'))
    db.update_value(5,'reactant1_ratio', 10)
    print(db.get_row_by_id(5))
    db.generate_and_store_product()
    print(db.get_row_by_id(241))
    
    properties_list = ['group_smiles', 'mol_functional_group']
    db.fill_data_properties(properties_list)
    
    db.export_to_csv('../data/database/reaction.csv')

    # nested_data = [(1, 2), (3, 4), (5, 6), (7, 8)]
    # encoded_str = encode_nested_structure(nested_data)
    # print(encoded_str)
    # decoded_data = decode_nested_structure(encoded_str)
    # print(decoded_data[0][1])
