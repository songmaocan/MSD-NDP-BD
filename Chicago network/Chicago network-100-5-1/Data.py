class Read_data:
    def __init__(self,multiplier):
        self.multiplier=multiplier

    def read_nodes(self):
        self.node_list=[]
        self.g_number_of_nodes=0
        with open("nodes.txt","r") as fl:
            lines=fl.readlines()
            for line in lines[1:]:
                str_list=line.strip().split("\t")
                node=Node()
                node.node_id=self.g_number_of_nodes
                self.node_list.append(node)
                self.g_number_of_nodes+=1

    def read_links(self):
        self.link_list = []
        self.g_number_of_links=0
        #1.read available links
        with open("links.txt","r") as fl:
            lines=fl.readlines()
            for line in lines[1:]:
                str_list = line.strip().split("\t")
                link=Link()
                link.link_id=self.g_number_of_links
                link.link_type=0
                link.from_node_id=int(str_list[1])-1
                link.to_node_id=int(str_list[2])-1
                link.travel_time_mean=float(str_list[3])
                link.travel_time_variance=(float(str_list[4]))**2
                self.link_list.append(link)
                self.g_number_of_links+=1

                self.node_list[link.from_node_id].outbound_nodes_list.append(link.to_node_id)
                self.node_list[link.from_node_id].outbound_links_list.append(link)
                self.node_list[link.from_node_id].outbound_nodes_number = len(self.node_list[link.from_node_id].outbound_nodes_list)

                self.node_list[link.to_node_id].inbound_nodes_list.append(link.from_node_id)
                self.node_list[link.to_node_id].inbound_links_list.append(link)
                self.node_list[link.to_node_id].inbound_nodes_number = len(self.node_list[link.to_node_id].inbound_nodes_list)


    def read_candidate_links(self):
        self.read_nodes()
        self.read_links()

        with open("OD pairs.txt","r") as fl:
            self.OD_pair_list=[]
            self.g_number_of_ODs=0
            lines=fl.readlines()
            for line in lines[1:]:
                str_list = line.strip().split("\t")
                od_pair=(int(str_list[1])-1,int(str_list[2])-1)
                self.OD_pair_list.append(od_pair)
                self.g_number_of_ODs+=1

        self.candidate_link_list=[]
        self.candidate_link_id_list = []
        self.g_number_of_candidate_links=0
        with open("Candidate links.txt","r") as fl:
            lines = fl.readlines()
            for line in lines[1:]:
                str_list = line.strip().split("\t")
                link = Link()
                link.link_id = self.g_number_of_links
                link.link_type = 1
                link.from_node_id = int(str_list[1]) - 1
                link.to_node_id = int(str_list[2]) - 1
                link.travel_time_mean = float(str_list[3])
                link.travel_time_variance = (float(str_list[4])) ** 2
                link.base_profit_for_lagrangian = [self.multiplier] * self.g_number_of_ODs
                self.link_list.append(link)
                self.candidate_link_list.append(link)
                self.candidate_link_id_list.append(link.link_id)
                self.g_number_of_links+=1
                self.g_number_of_candidate_links+=1

                self.node_list[link.from_node_id].outbound_nodes_list.append(link.to_node_id)
                self.node_list[link.from_node_id].outbound_links_list.append(link)
                self.node_list[link.from_node_id].outbound_nodes_number = len(self.node_list[link.from_node_id].outbound_nodes_list)

                self.node_list[link.to_node_id].inbound_nodes_list.append(link.from_node_id)
                self.node_list[link.to_node_id].inbound_links_list.append(link)
                self.node_list[link.to_node_id].inbound_nodes_number = len(self.node_list[link.to_node_id].inbound_nodes_list)

        return self.node_list,self.link_list,self.candidate_link_list,self.OD_pair_list,\
               self.g_number_of_links,self.g_number_of_nodes,self.g_number_of_candidate_links,self.g_number_of_ODs


class Node:
    def __init__(self):
        # self.node_index=None    #From 0
        self.node_id=None  #from 0
        self.outbound_nodes_list=[]
        self.outbound_nodes_number=0
        self.outbound_links_list=[]
        self.inbound_nodes_list = []
        self.inbound_nodes_number=0
        self.inbound_links_list = []


class Link:
    def __init__(self):
        self.link_id=None
        self.link_type=None  #1 candidate link ; 0 o.w.
        self.construction_Flag=1 #according to the KS
        self.from_node_id=None
        self.to_node_id=None
        self.travel_time_mean=None
        self.travel_time_variance=None
        self.base_profit_for_lagrangian = None

# class OD:
#     def __init__(self):
#         self.index=None
#         self.origin=None
#         self.destination=None