void swSPMV_host::SplitMeshMetis()
{
  int K = A.GetM();
  metis_int_t nb_part = 64;
  metis_int_t nb_elt = K, ncon = 1;
  metis_int_t total_comm;
  Seldon::Vector<metis_int_t> xadj(K + 1);
  xadj(0) = 0;
  for (int num_elm1 = 0; num_elm1 < K; num_elm1++)
  {
    xadj(num_elm1 + 1) = xadj(num_elm1);
    xadj(num_elm1 + 1) += (A.GetRowSize(num_elm1) - 1);  // exclude own element
  }

  Seldon::Vector<metis_int_t> adjncy(xadj(K));
  int nb = 0, num_elem2;
  for (int num_elm1 = 0; num_elm1 < K; num_elm1++)
  {
    for (int n = 0; n < A.GetRowSize(num_elm1); n++)
    {
      const int num_elem2 = A.GetIndex(num_elm1)[n];

      if (num_elem2 != num_elm1)
      {
        adjncy(nb) = num_elem2;
        nb++;
      }
    }
  }

  Seldon::Vector<metis_int_t> weight_element(K);
  weight_element.Fill(1.2);
  epart.Reallocate(K);

  METIS_PartGraphKway(&nb_elt,
                      &ncon,
                      xadj.GetData(),
                      adjncy.GetData(),
                      weight_element.GetData(),
                      NULL,
                      NULL,
                      &nb_part,
                      NULL,
                      NULL,
                      NULL,
                      &total_comm,
                      epart.GetData());

  Seldon::Vector<metis_int_t> col_permutation(K);
  col_permutation.Fill();
  Seldon::Vector<metis_int_t> row = col_permutation;
  Seldon::Vector<metis_int_t> epart2 = epart;
  Seldon::Sort(epart, col_permutation);
  Seldon::Matrix<double, Seldon::General, Seldon::ArrayRowSparse> B2(K, K);

  for (int num_elm1 = 0; num_elm1 < K; num_elm1++)
  {
    for (int n = 0; n < A.GetRowSize(num_elm1); n++)
    {
      const int num_elem2 = A.GetIndex(num_elm1)[n];
      const double coeff = A(num_elm1, num_elem2);
      if (coeff != 0) B2.AddInteraction(num_elm1, num_elem2, coeff);
    }
  }

  Seldon::ApplyPermutation(B2, row, col_permutation);
  Seldon::ApplyPermutation(B2, col_permutation, row);
  resA = B2;

  X_RowIndex2Thread = epart;

  Seldon::WriteMatrixMarket(resA, "matrixB3.mtx");

  Seldon::Vector<double> x(K);
  Seldon::Vector<double> b(K);
  x.Fill();
  Seldon::Mlt(A, x, b);
  Seldon::Vector<double> x2(K);
  Seldon::Vector<double> b2(K);
  x2.Fill();

  col_permutation.Fill();
  Seldon::Sort(epart2, col_permutation, x2);

  std::setprecision(15);

  double start_t = timestamp();
  Seldon::Mlt(B2, x2, b2);
  double end_t = timestamp();
  DISP(end_t - start_t);
  X = x2;
  resB = b2;
  cout << std::setprecision(20) << *max_element(b.GetData(), b.GetData() + K)
       << " vs " << *min_element(b2.GetData(), b2.GetData() + K) << endl;

  std::sort(b.GetData(), b.GetData() + K);
  std::sort(b2.GetData(), b2.GetData() + K);
  for (int i = 0; i < K; i++)
  {
    //        cout << (b(i) - b2(i)) / b(i);
    double error = std::abs((b(i) - b2(i)) / b(i));
    if (error > 1e-12)
    {
      DISP(error);
    }
  }
}
