




class Channel
{
public: 
	int NChannel;
	double Time[1024];
	double Volt[1024];
	
};


class Evento
{
public:
	int NEvento;
	std::vector<Channel> channel;
};

