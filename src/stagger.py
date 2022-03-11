import numpy as np 
from scipy.io import FortranFile

class stag(object):
    def __init__(self, folder='/mnt/beegfs/gemini/groups/bergemann/users/hoppe/Stagger/grid/t5777g44m00/', \
                 model='t5777g44m0005', snap=20):
        
        s = self
        
        s.mesh_file = folder + model + '.msh'
        s.file = folder + model + '_{:05}'.format(snap)
        s.data_file = s.file + '.dat'
        s.aux_file = s.file + '.aux'
        
        s.gz=5
        
        msh = FortranFile(s.mesh_file, 'r', header_dtype='>i4')
        s.nx, s.ny, s.nz = msh.read_record('>i4')
        s.nx = int(s.nx)
        s.ny = int(s.ny)
        s.nz = int(s.nz)
        
        _, _, s.xx, _, _, _ = np.split(msh.read_record('>f4'),6)
        _, _, s.zz, _, _, _ = np.split(msh.read_record('>f4'),6)
        _, _, s.yy, _, _, _ = np.split(msh.read_record('>f4'),6)
        s.zz = s.zz[s.gz:-s.gz]
        s.zz = -s.zz[::-1]

        msh.close()
        
        dat = FortranFile(s.data_file, 'r', header_dtype='>i4')
        _ = dat._read_size()
        datmem = np.memmap(dat._fp,'>f4', mode='r', order='F').reshape((-1, s.nx, s.nz, s.ny))
        datmem = np.transpose(datmem, axes=(0,1,3,2))
        datmem = datmem[..., s.gz:-s.gz]
        datmem = datmem[..., ::-1]
        s.rho  = datmem[0]
        s.ipx  = datmem[1]
        s.ipz  = datmem[2]
        s.ipy  = datmem[3]
        s.iie  = datmem[4]
        s.temp = datmem[5]
        dat.close()
        
        aux = FortranFile(s.aux_file, 'r', header_dtype='>i4')
        _ = aux._read_size()
        
        auxmem = np.memmap(aux._fp,'>f4', mode='r', order='F').reshape((-1, s.nx, s.nz, s.ny))
        auxmem = np.transpose(auxmem, axes=(0,1,3,2))
        auxmem = auxmem[..., s.gz:-s.gz]
        auxmem = auxmem[..., ::-1]
        s.lpp     = auxmem[0]
        s.lross   = auxmem[1]
        s.ltemp   = auxmem[2]
        s.lne     = auxmem[3]
        s.lplanck = auxmem[4]
        s.ltau    = auxmem[5]
        aux.close()
        
        self.nz = s.zz.size
        
    def get_vz(self, iz):
        c = 3./256.
        b = -25./256.
        a = 0.5-b-c
        
        ipz = self.ipz
        
        pz = a*(ipz[:,:,iz  ] + ipz[:,:,iz+1]) \
           + b*(ipz[:,:,iz-1] + ipz[:,:,iz+2]) \
           + c*(ipz[:,:,iz-2] + ipz[:,:,iz+3])
        
        vz = -pz / self.rho[:,:,iz]
        
        return vz

    def read_mem(self, field):
        return getattr(self,field)[:,:,:]