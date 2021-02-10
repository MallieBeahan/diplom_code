/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   moya.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/10/22 22:04:50 by Alexandr          #+#    #+#             */
/*   Updated: 2021/01/31 20:05:30 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */
#ifndef GLOBAL_VARS_H
# define GLOBAL_VARS_H

struct GlobalVars{
    double *Vx;
    double *Vy;
    double *Vz;
    double *coordx;
    double *coordy;
    double *coordz;
    double *Fx;
    double *Fy;
    double *Fz;
    double *Epot;
};

#endif