$(document).ready(function(){
    $(".owl-carousel").owlCarousel({
        loop: true,
        margin: 30,
        nav: true,
        singleItem: true,
        items: 1 // THIS IS IMPORTANT
        // responsive : {
        //     480:{
        //         items:1 // from zero to 480 width screen show 4 items
        //     },
        //     768:{
        //         items:3 // from zero to 480 width screen show 6 items
        //     },
        //     1024: {
        //         items: 6 // from zero to 480 width screen show 8 items
        //     }
        // }
    });
})