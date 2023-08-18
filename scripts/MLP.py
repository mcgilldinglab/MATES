import torch
import torch.nn.functional as F
from torch import nn

class MultiLayerPerceptron(nn.Module):
    def __init__(self, fam_num, input_dim, bias = True):
        super(MultiLayerPerceptron,self).__init__()
        self.n_input = input_dim
        self.n_fam = fam_num
        self.mlp = nn.Sequential(
            nn.Linear(self.n_input+self.n_fam , 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )
        
    def forward(self, embedding, BATCH_data, BATCH_SIZE):
        embedding = embedding.reshape(BATCH_SIZE,self.n_input)
        batch_id_encode = torch.eye(self.n_fam)[BATCH_data.type(torch.LongTensor)].view(BATCH_SIZE, self.n_fam).to('cuda:0')
        Y = torch.cat([embedding,batch_id_encode],axis = 1).reshape(BATCH_SIZE,1,embedding.shape[1]+self.n_fam)
        alpha = self.mlp(Y)
        return alpha

class MLP_loss(nn.Module):
    def __init__(self):
        super().__init__()
    def forward(self, alpha, TE_region, BATCH_SIZE, BIN_SIZE):
        alpha = alpha.reshape(BATCH_SIZE)
        loss = torch.zeros(BATCH_SIZE, dtype=torch.float)
        tmp = TE_region
        # [U_u_sum, U_m_sum, M_u_sum, M_m_sum, multi_region_num]
        #Multi = (tmp[:,-2] - tmp[:,-3]) * alpha
        #Unique = (tmp[:,1] - tmp[:,0])
        Multi = (tmp[:,-2] - tmp[:,-3]) * alpha / (tmp[:,-1]*BIN_SIZE)
        Unique = (tmp[:,1] - tmp[:,0]) / (tmp[:,-1]*BIN_SIZE)
        Multi = Multi.reshape(BATCH_SIZE)
        Unique = Unique.reshape(BATCH_SIZE)

        loss = torch.abs(Multi-Unique)
        return loss
    
