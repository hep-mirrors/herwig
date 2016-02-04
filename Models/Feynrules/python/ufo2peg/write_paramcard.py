# modified from the UFO version

class ParamCardWriter(object):
    
    def __init__(self, all_parameters):
        """write a valid param_card.dat"""
        
        list_of_parameters = [
              param 
              for param in all_parameters 
              if param.nature=='external'
            ]
        
        self.output = []
        
        self.write_card(list_of_parameters)
    
    def write_card(self, all_ext_param):
        # list all lhablock
        all_lhablock = set([param.lhablock for param in all_ext_param])
        
        # ordonate lhablock alphabeticaly
        list(all_lhablock).sort()
        
        for lhablock in all_lhablock:
            self.write_block(lhablock)
            [self.write_param(param, lhablock) for param in all_ext_param if \
                                                     param.lhablock == lhablock]
    def write_block(self, name):
        """ write a comment for a block"""
        if name!='DECAY':
            self.output.append('<< "Block %s\\n"' % name)

    def write_param(self, param, lhablock):
        
        lhacode=' '.join(['%3s' % key for key in param.lhacode])
        if lhablock != 'DECAY':
            text = '<< "  %s " << %s() << " # %s\\n"' % (lhacode, param.name, param.name ) 
        else:
            text = '<< "DECAY %s " << %s() << "\\n"' % (lhacode, param.name)
        self.output.append(text) 
