import torch

# Model size of Spatial-DC
def load_distribution_model(load_model_path):
    
    checkpoint = torch.load(load_model_path, map_location=torch.device('cpu'))
    return(checkpoint)
    # self.distriubtion_model.load_state_dict(checkpoint['model'])   


def cal_model_size(cpt):
    # 计算总大小（以字节为单位）
    total_size = 0
    state_dict = cpt["model"]
    for param_key, param_value in state_dict.items():
        # 获取参数的形状
        param_shape = param_value.size()
        # 计算参数的大小（以字节为单位）
        param_size = torch.numel(param_value) * param_value.element_size()
        # 累加大小
        total_size += param_size
    
    # 打印总大小（以MB为单位）
    total_size_MB = total_size / (1024 ** 2)
    print(f"Model size: {total_size_MB:.2f} MB")

load_model_path = "model_epoch200.pt"
cpt = load_distribution_model(load_model_path)
cal_model_size(cpt)
# Model size: 5.28 MB
