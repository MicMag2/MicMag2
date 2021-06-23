import magnum.module as module
import magnum.magneto as magneto



class Topology(module.Module):
    def __init__(self):
        super(Topology, self).__init__()
        setattr(self,"topo_method","berg-luescher")
        self.reset_state = False

    def calculates(self):
        return ["Q","Q_density"]

    def params(self):
        return ["topo_method"]

    def properties(self):
        return dict()

    def on_param_update(self,id):
        if id == "topo_method":
            self.reset_state = True
            #print('Q density should be resetted')

    def initialize(self, system):
        self.system = system

    def calculate(self, state, id):
        cache = state.cache
        if self.reset_state == True:
            if hasattr(cache,"Q"): del cache.Q
            if hasattr(cache,"Q_density"): del cache.Q_density
            self.reset_state = False
            #print("Reset Density")

        if id == "Q":
            if hasattr(cache,"Q"): return cache.Q
            if getattr(self,"topo_method") == "berg-luescher":
                cache.Q = magneto.topology_charge_berg_luescher(0,state.M,state.Ms) 
                return cache.Q
            elif getattr(self,"topo_method") == "berg-luescher-dual-lattice":
                cache.Q = magneto.topology_charge_berg_luescher_dual_lattice(0,state.M,state.Ms) 
                return cache.Q
            elif getattr(self,"topo_method") == "continuous":
                cache.Q = magneto.topology_charge_continuous(0,state.M,state.Ms)
                return cache.Q
            else:
                raise KeyError(f"Topology: topo_method = {getattr(self,'topo_method')} not defined")

        elif id == "Q_density":
            if hasattr(cache,"Q_density"): return cache.Q_density
            if getattr(self,"topo_method") == "berg-luescher":
                cache.Q_density = magneto.topology_charge_berg_luescher_density(state.M,state.Ms)
                return cache.Q_density
            elif getattr(self,"topo_method") == "berg-luescher-dual-lattice":
                cache.Q_density = magneto.topology_charge_berg_luescher_density_dual_lattice(state.M,state.Ms)
                return cache.Q_density
            elif getattr(self,"topo_method") ==  "continuous":  
                cache.Q_density = magneto.topology_charge_density_continuous(state.M,state.Ms)
                return cache.Q_density
            else:
                raise KeyError(f"Topology: topo_method = {getattr(self,'topo_method')} not defined")
       
        else: 
            raise KeyError(f"Topology: Can't calculate  {id}")
