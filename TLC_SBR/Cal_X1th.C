{
    double x1 = (280.0 + 290.0)/2.0;
    double x2 = (628.0 + 638.0)/2.0;

    double offset1 = (27.0 + 35.5)/2.0;
    double offset2 = (25.0 + 33.5)/2.0;
    
    double gain1 = 4.0/abs(27.0 - 35.5);
    double gain2 = 4.0/abs(25.0 - 33.5);
    
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
    
    std::cout << "X1thCal " << -c_offset << " " << -m_offset << " " << c_gain << " " << m_gain << std::endl;
}
