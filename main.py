import logging
import os
from src.database import DatabaseModule
from src.manager import run_parallel, run_simulation

logging.basicConfig(level=logging.INFO)
def main():
    db_PATH = 'data/database/reaction.db'
    try:
        if not os.path.exists(db_PATH):
            # Step 1: 建立数据库，并实例化
            db = DatabaseModule(db_path = db_PATH)

            # Step 2: 初始化数据库：填充反应条件、反应产物、分子特征等
            db.data_init()
            db.close()
        
        # Step 3: 并行运行数据库每行的反应任务获取结果
        # run_parallel()
        run_simulation(66)
    except Exception as e:
        logging.exception(f"启动失败: {e}")



if __name__ == "__main__":
    main()