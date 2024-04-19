import torch
import torch.nn.functional as F
from torch import nn

class AutoEncoder(nn.Module):
    def __init__(self,input_dim, output_dim, fam_num, bias = True):
        super(AutoEncoder,self).__init__()    
        self.n_input = input_dim
        self.n_output = output_dim
        self.n_fam = fam_num
        
        ##Encoder
        self.encoder = nn.Sequential(
            nn.Linear(self.n_input + self.n_fam, 1024, bias = bias),
            nn.ReLU(),
            nn.Linear(1024, 512, bias = bias),
            nn.ReLU(),
            nn.Linear(512, self.n_output, bias = bias),
            nn.ReLU(),
        )
        
        ##Decoder
        self.decoder = nn.Sequential(
            nn.Linear(self.n_output + self.n_fam, 512, bias = bias),
            nn.ReLU(),
            nn.Linear(512, 1024, bias = bias),
            nn.ReLU(),
            nn.Linear(1024, self.n_input, bias = bias),
            nn.ReLU(), 
        )
  
    def forward(self, TE_data, BATCH_data, BATCH_SIZE, ):
        reshaped_TE=torch.reshape(TE_data, (BATCH_SIZE,2001)).to('cuda:0')
        ##one-hot encoding TE Fam info
        batch_id_encode = torch.eye(self.n_fam)[BATCH_data.type(torch.LongTensor)].view(BATCH_SIZE, self.n_fam).to('cuda:0')
              
        new_reshaped_TE = torch.cat([reshaped_TE, batch_id_encode],axis = 1).reshape(BATCH_SIZE,1,self.n_input+self.n_fam)
        embedding = self.encoder(new_reshaped_TE).reshape(BATCH_SIZE,self.n_output)
        
        X = torch.cat([embedding,batch_id_encode],axis = 1).reshape(BATCH_SIZE,1,self.n_output+self.n_fam)
        
        reconstruction = self.decoder(X).reshape(BATCH_SIZE,self.n_input)
        return embedding, reconstruction
