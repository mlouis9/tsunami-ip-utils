from tsunami_ip_utils import RegionIntegratedSdfReader

model_1 = RegionIntegratedSdfReader('sphere_model_1.sdf')
model_1.convert_to_dict()
print(model_1.sdf_data['u-235'].keys())