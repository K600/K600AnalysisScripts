{
    double x1 = (261.0 + 271.0)/2.0;
    double x2 = (637.0 + 647.0)/2.0;

    double offset1 = (32.0 + 42.0)/2.0;
    double offset2 = (30.0 + 40.0)/2.0;
    
    double gain1 = 4.0/abs(32.0 - 42.0);
    double gain2 = 4.0/abs(30.0 - 40.0);
    
    double m_offset = (offset2-offset1)/(x2-x1);
    double c_offset = offset1 - (m_offset*x1);

    double m_gain = (gain2-gain1)/(x2-x1);
    double c_gain = gain1 - (m_gain*x1);
    
    /*
    std::cout << "m_offset: " << m_offset << std::endl;
    std::cout << "c_offset: " << c_offset << std::endl;

    std::cout << "m_gain: " << m_gain << std::endl;
    std::cout << "c_gain: " << c_gain << std::endl;
    */
    
    std::cout << "U1thCal " << -c_offset << " " << -m_offset << " " << c_gain << " " << m_gain << std::endl;
}
